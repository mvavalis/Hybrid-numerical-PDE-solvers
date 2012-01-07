#ifndef _LAPLACE_DD_H_
#define _LAPLACE_DD_H_

#include <pthread.h>
#include "dd.h"
#include "Problem.h"
#include "LaplaceSolve.h"
#include <MBA.h>
#include <UCButils.h>
//#include <base/point.h>
#include "more/euclid_norm.h"

/**
 * LaplaceDriver.h handles the last step of the domain decomposition.
 * LaplaceDriver<int dim> carries all the information needed to drive the
 * domain decomposition, and it also carries all the data needed by
 * the various -to be created- LaplaceSolve<int dim> objects in
 * order to solve their part of the problem.
 *
 * The way it goes:
 * LaplaceDriver<int dim> creates LaplaceSolve<int dim> objects and
 * gives them access to the required data, along with the information
 * needed to handle this data properly (through struct lapjob).
 */

//========================== 2D or 3D ==========================
template <int dim>
class LaplaceDriver
{
private:
	Problem<dim> &Prob;
	const struct pdd_prm<dim> &pp;
	
	static void* lapmain(void *arg);
public:
	LaplaceDriver(Problem<dim> &Prob, const struct pdd_prm<dim> &pp);
	void laplace_driver(double **node_coord, double *mc_sol_est);
	void laplace_driver(UCBspl::SplineSurface *dec_pl);
	
#if TEST_INTERFACE == 1
	void test_interfaces(double **node_coord, struct dec_hpl<2> *hpl);
#endif
};

template <int dim>
struct lapjob {
	int nof_jobs;
	int offset;
	struct lapprm<dim> *lapp;
};

template <int dim>
LaplaceDriver<dim>::LaplaceDriver(Problem<dim> &Prob, const struct pdd_prm<dim> &pp) : Prob(Prob), pp(pp) {}

template <int dim>
void* LaplaceDriver<dim>::lapmain(void *arg)
{
	struct lapjob<dim> *lapj = (struct lapjob<dim> *)arg;
	
	for (int i=0; i<lapj->nof_jobs; i++) {
		int job_id = lapj->offset + i;
		LaplaceSolve_main<dim>(lapj->lapp, job_id);
	}
	pthread_exit(NULL);
}

//============================= 2D =============================
bool solution_required(const struct solve_prm *sp, int xi, int yi)
{
	if (!sp->selective_sol) {return true;}
	
	int i;
	for (i=0; i<sp->nof_subd; i++) {
		if (sp->subp[i].subd_coord[Xd] == xi && sp->subp[i].subd_coord[Yd] == yi) {
			return true;
		}
	}
	return false;
}

/* create INTER2D_NOF_TEST_NODES nodes in \omega and compare to the real solution */
#define INTER2D_NOF_TEST_NODES 2000
template <>
void LaplaceDriver<2>::test_interfaces(double **node_coord, struct dec_hpl<2> *hpl)
{
	int Xlines = pp.dp.Xlines;
	int Ylines = pp.dp.Ylines;
	double *Ycoord = pp.dp.Ycoord,  *Xcoord = pp.dp.Xcoord;
	
	double test_node_coord[INTER2D_NOF_TEST_NODES][2];
	double test_est_values[INTER2D_NOF_TEST_NODES];
	double test_actual_values[INTER2D_NOF_TEST_NODES];
	
	if (!(Ylines > 0 || Xlines > 0)) {
		return;
	}
	
	struct norm norm;
	std::cout<<std::endl<<"*========================================*"<<std::endl;
	std::cout<<           "* Interface test -- nof test nodes: "<<INTER2D_NOF_TEST_NODES<<" *"<<std::endl;
	std::cout<<           "*========================================*"<<std::endl;
	for (int j=0; j<Ylines+Xlines; j++) {
		for (int i=0; i<INTER2D_NOF_TEST_NODES; i++) {
	// 		test_node_coord[i] = new double[2];
			if (j < Ylines) {
				test_node_coord[i][Xd] = Ycoord[j];
				test_node_coord[i][Yd] = pp.D[Yd] * (i /(double) INTER2D_NOF_TEST_NODES);
			}
			else {
				test_node_coord[i][Xd] = pp.D[Xd] * (i /(double) INTER2D_NOF_TEST_NODES);
				test_node_coord[i][Yd] = Xcoord[j-Xlines];
			}
			
#if   INTERP2D_LIB == 1
			test_est_values[i] =
			    spline_b_val(hpl[j].nof_nodes,  hpl[j].spl_node_coord,  hpl[j].mc_sol_est,
			                 test_node_coord[i][(j < Ylines) ? Yd:Xd]);
#elif INTERP2D_LIB == 2
			test_est_values[i] = gsl_spline_eval(hpl[j].spline, test_node_coord[i][(j < Ylines) ? Yd:Xd], hpl[j].acc);
#endif
			
			test_actual_values[i] = Prob.test_u(test_node_coord[i]);
		}
		
		std::cout<<"Interface #"<<j<<" under test (nof interface nodes: "<<hpl[j].nof_nodes<<"):"<<std::endl;
// 		for (int i=0; i<hpl[j].nof_nodes; i++) {
// 			std::cout<<hpl[j].spl_node_coord[i]<<" ";
// 			std::cout<<": "<<hpl[j].mc_sol_est[i]<<std::endl;
// 		}
// 		std::cout<<"-----------"<<std::endl;
		
		norm = calc_norm(test_est_values, test_actual_values, INTER2D_NOF_TEST_NODES);
		cout << "Euclidean norm: " << norm.diff_norm << endl;
		cout << "Relative euclidean norm: " << norm.diff_norm/norm.actual_norm << endl;
		std::cout<<"-----------"<<std::endl;
	}
}

//Make ready struct lapjob<2> s and create pthreads (pthread_create(,, lapmain, lapj))
template <>
void LaplaceDriver<2>::laplace_driver(double **node_coord, double *mc_sol_est)
{
	int i, j, nof_threads, nof_jobs;
	int nof_jpthread, ntjp; //nof jobs per thread, nof threads with plus one job
	pthread_t *pmc;
	struct lapjob<2> *lapj;
	struct dec_hpl<2> *hpl;
	
	
	int Xlines = pp.dp.Xlines;  int nodes_X = pp.dp.nodes_X;
	int Ylines = pp.dp.Ylines;  int nodes_Y = pp.dp.nodes_Y;
	double *Ycoord = pp.dp.Ycoord,  *Xcoord = pp.dp.Xcoord;
	
	nof_jobs = (pp.sp.selective_sol == false) ? (Xlines+1)*(Ylines+1) : pp.sp.nof_subd;
	
	nof_threads = (pp.max_nof_threads > nof_jobs) ? nof_jobs:pp.max_nof_threads;
	nof_jpthread = nof_jobs/nof_threads;
	ntjp = nof_jobs - nof_jpthread*nof_threads;
	
	std::cout<<"\n*** Laplace ***"<<std::endl
			<<"nof threads:\t"<<nof_threads<<std::endl
			<<"nof jobs:\t"<<nof_jobs<<std::endl;
	
	
	double *Ylines_coordY = new double[nodes_Y];
	double *Xlines_coordX = new double[nodes_X];
	int allnodes_Y = Ylines*nodes_Y;
	for (i=0; i<nodes_Y; i++) {Ylines_coordY[i] = node_coord[i][Yd];}
	for (i=0; i<nodes_X; i++) {Xlines_coordX[i] = node_coord[allnodes_Y + i][Xd];}
	
	hpl = new struct dec_hpl<2>[Xlines + Ylines];
	for (i=0; i<Ylines; i++) {
		hpl[i].nof_nodes = nodes_Y;
#if   INTERP2D_LIB == 1
		hpl[i].spl_node_coord = Ylines_coordY;
		hpl[i].mc_sol_est     = &mc_sol_est[i*nodes_Y];
#elif INTERP2D_LIB == 2
		hpl[i].acc    = gsl_interp_accel_alloc();
		hpl[i].spline = gsl_spline_alloc(gsl_interp_cspline, nodes_Y);
		gsl_spline_init(hpl[i].spline, Ylines_coordY, &mc_sol_est[i*nodes_Y], nodes_Y);
#endif
	}
	for (i=0; i<Xlines; i++) {
		hpl[Ylines + i].nof_nodes = nodes_X;
#if   INTERP2D_LIB == 1
		hpl[Ylines + i].spl_node_coord = Xlines_coordX;
		hpl[Ylines + i].mc_sol_est     = &mc_sol_est[allnodes_Y + i*nodes_X];
#elif INTERP2D_LIB == 2
		hpl[Ylines + i].acc    = gsl_interp_accel_alloc();
		hpl[Ylines + i].spline = gsl_spline_alloc(gsl_interp_cspline, nodes_X);
		gsl_spline_init(hpl[Ylines + i].spline, Xlines_coordX, &mc_sol_est[allnodes_Y + i*nodes_X], nodes_X);
#endif
	}
// #define INTERFACE_TEST 0
// #if INTERFACE_TEST == 1 /* test the 1st boundary */
// 		//#include <stdlib.h>
// 		//#include <time.h>
// 		#define INTER2D_NOF_TEST_NODES 1000
// 		#define SPL1 0
// 		/* create INTER2D_NOF_TEST_NODES nodes in \omega and compare to the real solution */
// 		double *test_node_coord[INTER2D_NOF_TEST_NODES];
// 		double test_est_values[INTER2D_NOF_TEST_NODES];
// 		//srand(time(NULL));
// 		if (Ylines > 0 || Xlines > 0) {
// 			for (int i=0; i<INTER2D_NOF_TEST_NODES; i++) {
// 				test_node_coord[i] = new double[2];
// 				if (Ylines > 0) {
// 					test_node_coord[i][Xd] = node_coord[0][Xd];
// 					test_node_coord[i][Yd] = pp.D[Yd] * /*(rand() /(double) RAND_MAX)*/ (i /(double) INTER2D_NOF_TEST_NODES);
// 				}
// 				else {
// 					test_node_coord[i][Xd] = pp.D[Xd] * /*(rand() /(double) RAND_MAX)*/ (i /(double) INTER2D_NOF_TEST_NODES);
// 					test_node_coord[i][Yd] = node_coord[0][Yd];
// 				}
// 				test_est_values[i] =
// 					spline_b_val(hpl[SPL1].nof_nodes, hpl[SPL1].spl_node_coord, hpl[SPL1].mc_sol_est, test_node_coord[i][(Ylines > 0) ? Yd:Xd]);
// 			}
// 			
// 			std::cout<<"1st interface under test (nof_nodes: "<<hpl[SPL1].nof_nodes<<"):"<<std::endl;
// 			for (int i=0; i<hpl[SPL1].nof_nodes; i++) {
// 				std::cout<<hpl[SPL1].spl_node_coord[i]<<" ";
// 				std::cout<<": "<<hpl[SPL1].mc_sol_est[i]<<std::endl;
// 			}
// 			std::cout<<"-----------"<<std::endl;
// 			
// 			std::cout<<"nof test nodes: "<<INTER2D_NOF_TEST_NODES<<std::endl;
// 			calc_norm<2>(test_node_coord, test_est_values, pp.D, INTER2D_NOF_TEST_NODES);
// 			std::cout<<"-----------"<<std::endl<<std::endl;
// 		}
// #endif
	
	std::vector<Point<2> > grid_hrect_p1, grid_hrect_p2; //start and end of hyper-rectangle
	for (i=0; i<Xlines+1; i++) { //y
		for (j=0; j<Ylines+1; j++) { //x
			if (solution_required(&pp.sp, j+1, i+1)) {
				double p1_x = (j == 0) ? 0. : Ycoord[j-1];
				double p1_y = (i == 0) ? 0. : Xcoord[i-1];
				double p2_x = (j == Ylines) ? pp.D[Xd] : Ycoord[j];
				double p2_y = (i == Xlines) ? pp.D[Yd] : Xcoord[i];
				cout<<"++++++"<<endl;
				cout<<"p1: "<<p1_x;
				cout<<" "<<p1_y<<endl;
				cout<<"p2: "<<p2_x;
				cout<<" "<<p2_y<<endl;
				cout<<"++++++"<<endl;
				grid_hrect_p1.push_back(Point<2>(p1_x, p1_y));
				grid_hrect_p2.push_back(Point<2>(p2_x, p2_y));
			}
		}
	}
	
	
	bool ***init_bound = new bool**[nof_jobs];
	int ***hp_spline   = new int**[nof_jobs];
	for (i=0; i<nof_jobs; i++) {
		init_bound[i] = new bool*[2/*dim*/];
		hp_spline [i] = new int*[2/*dim*/];
		for (j=0; j<2/*dim*/; j++) {
			init_bound[i][j] = new bool[2];
			hp_spline [i][j] = new int[2];
		}
	}
	
	/*
	 * Ylines are enumerated before Xlines 'cause the X dimension is
	 * traversed first.
	 * (./main.out 2D <sth> <sth>  4 3  <sth> <sth>  <sth> <sth>  <sth>)
	 *
	 *                  Y ^
	 *                    |
	 *                    |
	 *                    |----|----|----|----|
	 *                    |  8 |  9 | 10 | 11 |
	 * Xline (id: 4) ---> |----|----|----|----|
	 *                    |  4 |  5 |  6 |  7 |
	 * Xline (id: 3) ---> |----|----|----|----|
	 *                    |  0 |  1 |  2 |  3 |
	 *                    |----|----|----|----|-----> X
	 *                         ^    ^    ^
	 *                         |    |    |
	 * Yline (id: 0) ----------|    |    |
	 * Yline (id: 1) ---------------|    |
	 * Yline (id: 2) --------------------|
	 */
	
	int job_id = 0;
	for (i=0; i<Xlines+1; i++) { //y
		for (j=0; j<Ylines+1; j++) { //x
			if (solution_required(&pp.sp, j+1, i+1)) {
				//array indices: [job_id][dim][left or right]
				
				init_bound[job_id][Xd][0] = (j==0) ? true:false;
				init_bound[job_id][Xd][1] = (j==Ylines) ? true:false;
				init_bound[job_id][Yd][0] = (i==0) ? true:false;
				init_bound[job_id][Yd][1] = (i==Xlines) ? true:false;
				
				hp_spline [job_id][Xd][0] = (j-1);
				hp_spline [job_id][Xd][1] = j;
				hp_spline [job_id][Yd][0] = Ylines + (i-1);
				hp_spline [job_id][Yd][1] = Ylines + i;
				
				job_id++;
				
	/*
	// 	std::cout<<std::endl;
	// 	std::cout << "part: " << (job_id) << " (" << i << ", " << j << ")" << std::endl;
	// #define some_print(dim, right) \
	// do { \
	// 	if (init_bound[job_id][dim][right]) {std::cout << "-";} \
	// 	else {std::cout << hp_spline[job_id][dim][right];} \
	// 	std::cout << std::endl; \
	// } while (0)
	// 	std::cout << "Xd left : ";  some_print(Xd, 0);
	// 	std::cout << "Xd right: ";  some_print(Xd, 1);
	// 	std::cout << "Yd left : ";  some_print(Yd, 0);
	// 	std::cout << "Yd right: ";  some_print(Yd, 1);
	// #undef some_print*/
			}
		}
	}
	
	struct lapprm<2> lapp = {Prob, hpl, grid_hrect_p1, grid_hrect_p2, init_bound, hp_spline, &pp.sp};
	
	/*
// 	bool left; //left (or right?) bound
// 	int min_dim; // Xd, Yd, Zd
// 	double min;
// int job_id = 0;
// int dim = 2;
// const Point<2> p(1., 0.7);
// 	
// 	left = true;
// 	min = p(0)-lapp.grid_hrect_p1[job_id][0];  min_dim = 0;  left = true;
// 	for (i=1; i<dim; i++) {
// 		if (p(i)-lapp.grid_hrect_p1[job_id][i] < min) {min = p(i)-lapp.grid_hrect_p1[job_id][i];  min_dim = i;}
// 	}
// 	for (i=0; i<dim; i++) {
// 		if (lapp.grid_hrect_p2[job_id][i]-p(i) < min) {min = lapp.grid_hrect_p2[job_id][i]-p(i);  min_dim = i;  left = false;}
// 	}
// 	
// 	
// 	hpl = &lapp.hpl[lapp.hp_spline[job_id][min_dim][left ? 0:1]];
// 	std::cout << "job_id: "<< job_id <<std::endl;
// 	std::cout << "min_dim: "<< min_dim <<std::endl;
// 	std::cout << "left ? 0:1: "<< (left ? 0:1) <<std::endl;
// 	std::cout << "## lapp.hp_spline[job_id][min_dim][left ? 0:1]: "<< lapp.hp_spline[job_id][min_dim][left ? 0:1] <<std::endl;
// 	
// 	std::cout << "#### spline_b_val: "<<std::endl;
// 	for (i=0; i<hpl->nof_nodes; i++) {
// 		std::cout << hpl->spl_node_coord[i] << ": " << hpl->mc_sol_est[i] << std::endl;
// 	}
// 	std::cout << "#### spline_gg_val: "<<std::endl;
// 	std::cout << spline_b_val(hpl->nof_nodes, hpl->spl_node_coord, hpl->mc_sol_est, (min_dim==Xd) ? p[1]:p[0])
// 		<<std::endl;
// 	exit(0);
	*/
	
	
	pmc = new pthread_t[nof_threads];
	lapj = new lapjob<2>[nof_threads];
	
	std::cout<<"ntjp: "<<ntjp<<", nof_threads: "<<nof_threads<<std::endl;
	for (i=0   ; i<ntjp       ; i++) {
		lapj[i].nof_jobs = nof_jpthread+1;
		lapj[i].offset = i*(nof_jpthread+1);
		lapj[i].lapp = &lapp;
	}
	for (i=ntjp; i<nof_threads; i++) {
		lapj[i].nof_jobs = nof_jpthread;
		lapj[i].offset = ntjp*(nof_jpthread+1) + (i-ntjp)*nof_jpthread;
		lapj[i].lapp = &lapp;
	}
	
#if TEST_INTERFACE == 1
	test_interfaces(node_coord, hpl);
#endif
	
	for (i=0; i<nof_threads; i++) {pthread_create(&pmc[i], NULL, lapmain, &lapj[i]);}
	for (i=0; i<nof_threads; i++) {pthread_join(pmc[i], NULL);}
	
// 	for (i=0; i<nof_jobs; i++) {delete &(*laplace_solve[i]);}
// 	delete [] pmc;  delete [] lapj;
}

//============================= 3D =============================
bool solution_required(const struct solve_prm *sp, int xi, int yi, int zi)
{
	if (!sp->selective_sol) {return true;}
	
	int i;
	for (i=0; i<sp->nof_subd; i++) {
		if (sp->subp[i].subd_coord[Xd] == xi &&
			sp->subp[i].subd_coord[Yd] == yi &&
			sp->subp[i].subd_coord[Zd] == zi) {
			return true;
		}
	}
	return false;
}

template <>
void LaplaceDriver<3>::laplace_driver(UCBspl::SplineSurface *dec_pl)
{
	int i, j, k, nof_threads, nof_jobs;
	int nof_jpthread, ntjp; //nof jobs per thread, nof threads with plus one job
	pthread_t *pmc;
	struct lapjob<3> *lapj;
	struct dec_hpl<3> *hpl;
	bool ***init_bound;
	int ***hp_spline;
	
	
	int XYplanes = pp.dp.XYplanes[0];  int YZplanes = pp.dp.YZplanes[0];  int XZplanes = pp.dp.XZplanes[0];
	double *YZcoord = pp.dp.YZcoord,  *XZcoord = pp.dp.XZcoord,  *XYcoord = pp.dp.XYcoord;
	
	nof_jobs = (pp.sp.selective_sol == false) ? (XYplanes+1)*(YZplanes+1)*(XZplanes+1) : pp.sp.nof_subd;
	
	nof_threads = (pp.max_nof_threads > nof_jobs) ? nof_jobs:pp.max_nof_threads;
	nof_jpthread = nof_jobs/nof_threads;
	ntjp = nof_jobs - nof_jpthread*nof_threads;
	
	std::cout<<"\n*** Laplace ***"<<std::endl
			<<"nof threads:\t"<<nof_threads<<std::endl
			<<"nof jobs:\t"<<nof_jobs<<std::endl;
	
	hpl = new struct dec_hpl<3>[YZplanes + XZplanes + XYplanes];
	for (i=0; i<YZplanes+XZplanes+XYplanes; i++) {hpl[i].dec_pl = &dec_pl[i];}
	
	
	std::vector<Point<3> > grid_hrect_p1, grid_hrect_p2;
	for (i=0; i<XYplanes+1; i++) { //z
		for (j=0; j<XZplanes+1; j++) { //y
			for (k=0; k<YZplanes+1; k++) { //x
				if (solution_required(&pp.sp, k+1, j+1, i+1)) {
					double p1_x = (k == 0) ? 0. : YZcoord[k-1];
					double p1_y = (j == 0) ? 0. : XZcoord[j-1];
					double p1_z = (i == 0) ? 0. : XYcoord[i-1];
					double p2_x = (k == YZplanes) ? pp.D[Xd] : YZcoord[k];
					double p2_y = (j == XZplanes) ? pp.D[Yd] : XZcoord[j];
					double p2_z = (i == XYplanes) ? pp.D[Zd] : XYcoord[i];
					cout<<"++++++"<<endl;
					cout<<"p1: "<<p1_x;
					cout<<" "<<p1_y;
					cout<<" "<<p1_z<<endl;
					cout<<"p2: "<<p2_x;
					cout<<" "<<p2_y;
					cout<<" "<<p2_z<<endl;
					cout<<"++++++"<<endl;
					grid_hrect_p1.push_back(Point<3>(p1_x, p1_y, p1_z));
					grid_hrect_p2.push_back(Point<3>(p2_x, p2_y, p2_z));
				}
			}
		}
	}
	
	
	init_bound = new bool**[nof_jobs];
	hp_spline  = new int**[nof_jobs];
	for (i=0; i<nof_jobs; i++) {
		init_bound[i] = new bool*[3/*dim*/];
		hp_spline[i]  = new int*[3/*dim*/];
		for (j=0; j<3/*dim*/; j++) {
			init_bound[i][j] = new bool[2];
			hp_spline[i][j]  = new int[2];
		}
	}
	
	//determine init_bound[3/*dim*/][2] s and hp_spline[3/*dim*/][2] s
	int job_id = 0;
	for (i=0; i<XYplanes+1; i++) { //z
		for (j=0; j<XZplanes+1; j++) { //y
			for (k=0; k<YZplanes+1; k++) { //x
				if (solution_required(&pp.sp, k+1, j+1, i+1)) {
					//array indices: [job_id][dim][left or right]
					
					init_bound[job_id][Xd][0] = (k==0) ? true:false;
					init_bound[job_id][Xd][1] = (k==YZplanes) ? true:false;
					init_bound[job_id][Yd][0] = (j==0) ? true:false;
					init_bound[job_id][Yd][1] = (j==XZplanes) ? true:false;
					init_bound[job_id][Zd][0] = (i==0) ? true:false;
					init_bound[job_id][Zd][1] = (i==XYplanes) ? true:false;
					
					hp_spline [job_id][Xd][0] = (k-1);
					hp_spline [job_id][Xd][1] = k;
					hp_spline [job_id][Yd][0] = YZplanes + (j-1);
					hp_spline [job_id][Yd][1] = YZplanes + j;
					hp_spline [job_id][Zd][0] = YZplanes + XZplanes + (i-1);
					hp_spline [job_id][Zd][1] = YZplanes + XZplanes + i;
					
					job_id++;
					
	// 	std::cout<<std::endl;
	// 	std::cout << "part: " << (job_id) << " (" << i << ", " << j << ", " << k << ")" << std::endl;
	// #define some_print(dim, right)
	// do {
	// 	if (init_bound[job_id][dim][right]) {std::cout << "-";}
	// 	else {std::cout << hp_spline[job_id][dim][right];}
	// 	std::cout << std::endl;
	// } while (0)
	// 	std::cout << "Xd left : ";  some_print(Xd, 0);
	// 	std::cout << "Xd right: ";  some_print(Xd, 1);
	// 	std::cout << "Yd left : ";  some_print(Yd, 0);
	// 	std::cout << "Yd right: ";  some_print(Yd, 1);
	// 	std::cout << "Zd left : ";  some_print(Zd, 0);
	// 	std::cout << "Zd right: ";  some_print(Zd, 1);
	// #undef some_print
				}
			}
		}
	}
	
	
	struct lapprm<3> lapp = {Prob, hpl, grid_hrect_p1, grid_hrect_p2, init_bound, hp_spline, &pp.sp};
	
	pmc = new pthread_t[nof_threads];
	lapj = new lapjob<3>[nof_threads];
	
	std::cout<<"ntjp: "<<ntjp<<", nof_threads: "<<nof_threads<<std::endl;
	for (i=0   ; i<ntjp       ; i++) {
		lapj[i].nof_jobs = nof_jpthread+1;
		lapj[i].offset = i*(nof_jpthread+1);
		lapj[i].lapp = &lapp;
	}
	for (i=ntjp; i<nof_threads; i++) {
		lapj[i].nof_jobs = nof_jpthread;
		lapj[i].offset = ntjp*(nof_jpthread+1) + (i-ntjp)*nof_jpthread;
		lapj[i].lapp = &lapp;
	}
	
	for (i=0; i<(int)grid_hrect_p1.size(); i++) {
		std::cout<<"hyper_rect: "<<grid_hrect_p1[i](0)<<" "<<grid_hrect_p1[i](1)<<" "<<grid_hrect_p1[i](2)<<std::endl;
	}
	
	for (int i=0; i<nof_threads; i++) {pthread_create(&pmc[i], NULL, lapmain, &lapj[i]);}
	
	for (i=0; i<nof_threads; i++) {pthread_join(pmc[i], NULL);}
}

#endif
