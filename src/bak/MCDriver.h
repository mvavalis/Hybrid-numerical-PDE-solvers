#ifndef _MC_H_
#define _MC_H_

#include <iostream>
#include <cmath>
#include <ctime>
#include <pthread.h>
#include <gsl/gsl_rng.h>
#include "dd.h"
#include "more/Vdcbin.h"

template <int dim>
class MCDriver
{
private:
	struct output {
		//double mnof_steps; //mean of the nof_steps required reach the boundary
		double msol_est; //mean of the solution estimate
	};
	
	Problem<dim> &Prob;
	const struct pdd_prm<dim> &pp;
	
	gsl_rng *rng; //random number generator
	Vdcbin vdcX, vdcY; //van der Corput sequence
	double btol; //boundary tolerance
	double D[dim], btol_D[dim];  //\Omega: dimension length, and dimension length minus boundary tolerance
	
	double calc_sphere_rad(const double *x);
	double a(double d);
	void rand_update_x(double *x, double d);
	void quasi_update_x(double *x, double d);
	void rand_update_y(const double *x, double *y, double d);
	void quasi_update_y(const double *x, double *y, double d);
	
	void init_quasi(int size);
	double quasirand_uniform();
	
	struct output solve(int nof_walks, const double *x_start);
	static void* mcmain(void *arg);
public:
	MCDriver(Problem<dim> &Prob, const struct pdd_prm<dim> &pp);
	void monte_carlo(int nof_nodes, double** node_coord, double* &sol_est);
};

template <int dim>
struct mcjob {
	MCDriver<dim> *mc;
	int nof_walks;
	int nof_nodes;
	double **node_coord;
	double *sol_est;
};

//#include "more/Ranq1.h"
//Ranq1 ranq(time(NULL)); //random number generator
//#define rand_uniform() (ranq.doub())
#define rand_uniform() (gsl_rng_uniform(rng))

//============================= 2D =============================
template <> inline double MCDriver<2>::calc_sphere_rad(const double *x)
{
	double min;
	min = (D[0]-x[0]);
	if (x[0] < min) {min = x[0];}
	if (D[1]-x[1] < min) {min = D[1]-x[1];}
	if (x[1] < min) {min = x[1];}
	return min;
}

template <> inline double MCDriver<2>::a(double d)
{
	return d*d/4.;
}

template <> inline void MCDriver<2>::rand_update_x(double *x, double d)
{
	double angle = rand_uniform()*(2.*M_PI);
	x[0] += d*cos(angle);
	x[1] += d*sin(angle);
}

template <> inline void MCDriver<2>::quasi_update_x(double *x, double d)
{
	double angle = quasirand_uniform()*(2.*M_PI);
	x[0] += d*cos(angle);
	x[1] += d*sin(angle);
}

// template <> inline void MCDriver<2>::quasi_update_x(double *x, double d)
// {
// 	double angle = vdcX.doub_cz()*(2.*M_PI);
// 	x[0] += d*cos(angle);
// 	x[1] += d*sin(angle);
// }

template <> inline void MCDriver<2>::rand_update_y(const double *x, double *y, double d)
{
	double angle = rand_uniform()*(2.*M_PI);
	double rad;
	do {rad = rand_uniform()*d;} while (((4.*rad)/(d*d))*log(d/rad) < rand_uniform()*(4./(M_E*d)));
	y[0] = x[0] + rad*cos(angle);
	y[1] = x[1] + rad*sin(angle);
}

template <> inline void MCDriver<2>::quasi_update_y(const double *x, double *y, double d)
{
	double angle = quasirand_uniform()*(2.*M_PI);
	double rad;
	do {rad = quasirand_uniform()*d;} while (((4.*rad)/(d*d))*log(d/rad) < quasirand_uniform()*(4./(M_E*d)));
	y[0] = x[0] + rad*cos(angle);
	y[1] = x[1] + rad*sin(angle);
}

// template <> inline void MCDriver<2>::quasi_update_y(const double *x, double *y, double d)
// {
// 	double angle = vdcY.doub_cz()*(2.*M_PI);
// 	double rad;
// 	do {rad = rand_uniform()*d;} while (((4.*rad)/(d*d))*log(d/rad) < rand_uniform()*(4./(M_E*d)));
// 	y[0] = x[0] + rad*cos(angle);
// 	y[1] = x[1] + rad*sin(angle);
// }

//============================= 3D =============================
template <> inline double MCDriver<3>::calc_sphere_rad(const double *x)
{
	double min = (D[0]-x[0]);
	if (x[0] < min) {min = x[0];}
//  		double ppD[3] = {2., 2., 2.};
//  		Problem<3> Probb(ppD);
	if (D[1]-x[1] < min) {min = D[1]-x[1];}
	if (x[1] < min) {min = x[1];}
	if (D[2]-x[2] < min) {min = D[2]-x[2];}
	if (x[2] < min) {min = x[2];}
	return min;
}

template <> inline double MCDriver<3>::a(double d)
{
	return d*d/6.;
}

template <> inline void MCDriver<3>::rand_update_x(double *x, double d)
{
	double theta = rand_uniform()*(2.*M_PI);
	double phi = rand_uniform()*(M_PI);
	x[0] += d*sin(phi)*cos(theta);
	x[1] += d*sin(phi)*sin(theta);
	x[2] += d*cos(phi);
}

template <> inline void MCDriver<3>::quasi_update_x(double *x, double d)
{
	double theta = quasirand_uniform()*(2.*M_PI);
	double phi = quasirand_uniform()*(M_PI);
	x[0] += d*sin(phi)*cos(theta);
	x[1] += d*sin(phi)*sin(theta);
	x[2] += d*cos(phi);
}

// template <> inline void MCDriver<3>::quasi_update_x(double *x, double d)
// {
// 	double theta = vdcX.doub_cz()*(2.*M_PI);
// 	double phi = rand_uniform()*(M_PI);
// 	x[0] += d*sin(phi)*cos(theta);
// 	x[1] += d*sin(phi)*sin(theta);
// 	x[2] += d*cos(phi);
// }

template <> inline void MCDriver<3>::rand_update_y(const double *x, double *y, double d)
{
	double theta = rand_uniform()*(2.*M_PI);
	double rad, phi;
	do {rad = rand_uniform()*d;  phi = rand_uniform()*M_PI;}
	while ((3./(d*d*d))*rad*(d-rad)*sin(phi) < rand_uniform()*((3./4.)*d));
	y[0] = x[0] + rad*sin(phi)*cos(theta);
	y[1] = x[1] + rad*sin(phi)*sin(theta);
	y[2] = x[2] + rad*cos(phi);
}

template <> inline void MCDriver<3>::quasi_update_y(const double *x, double *y, double d)
{
	double theta = quasirand_uniform()*(2.*M_PI);
	double rad, phi;
	do {rad = quasirand_uniform()*d;  phi = quasirand_uniform()*M_PI;}
	while ((3./(d*d*d))*rad*(d-rad)*sin(phi) < quasirand_uniform()*((3./4.)*d));
	y[0] = x[0] + rad*sin(phi)*cos(theta);
	y[1] = x[1] + rad*sin(phi)*sin(theta);
	y[2] = x[2] + rad*cos(phi);
}

// template <> inline void MCDriver<3>::quasi_update_y(const double *x, double *y, double d)
// {
// 	double theta = vdcY.doub_cz()*(2.*M_PI);
// 	double rad, phi;
// 	do {rad = rand_uniform()*d;  phi = rand_uniform()*M_PI;}
// 	while ((3./(d*d*d))*rad*(d-rad)*sin(phi) < rand_uniform()*((3./4.)*d));
// 	y[0] = x[0] + rad*sin(phi)*cos(theta);
// 	y[1] = x[1] + rad*sin(phi)*sin(theta);
// 	y[2] = x[2] + rad*cos(phi);
// }


//========================== 2D or 3D ==========================
template <int dim>
MCDriver<dim>::MCDriver(Problem<dim> &Prob, const struct pdd_prm<dim> &pp) : Prob(Prob), pp(pp), vdcX(Vdcbin()), vdcY(Vdcbin())
{
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(rng, time(NULL));
	
	int i;
	for (i=0; i<dim; i++) {D[i] = Prob.D[i];}
	btol = pp.btol;
	for (i=0; i<dim; i++) {btol_D[i] = D[i]-btol;}
}

#if RAND_MODE == 2
static double *qarray;
static int qindex;
static int qsize;
template <int dim>
void MCDriver<dim>::init_quasi(int size)
{
	int tmp_index;
	double tmp;
	
	qsize = size;
	qindex = 0;
	qarray = new double[size];
	
	for (int i=0; i<size; i++) {
		qarray[i] = vdcX.doub_cz();
// 		qarray[i] = rand_uniform(); /*used for debugging*/
	}
	
// #define NUMMMM 107
// 	for (int i=0; i<size; i++) {
// 		for (int j=0; j<NUMMMM; j++) {
// 			tmp = qarray[(i+j)];
// 			tmp_index = (i+j) + gsl_rng_uniform_int(rng, (i+NUMMMM)-(i+j));
// 			
// 			qarray[(i+j)] = qarray[tmp_index];
// 			qarray[tmp_index] = tmp;
// 		}
// 		i+=NUMMMM-1;
// 	}
	
	for (int i=0; i<size; i++) {
		tmp = qarray[i];
		tmp_index = i + gsl_rng_uniform_int(rng, size-i);
		
		qarray[i] = qarray[tmp_index];
		qarray[tmp_index] = tmp;
	}
}

template <int dim>
double MCDriver<dim>::quasirand_uniform()
{
	if (qindex == qsize) {
		delete [] qarray;
		init_quasi(qsize);
		std::cout<<"init_quasi(qsize);"<<std::endl;
	}
	
	return qarray[qindex++];
}
#endif

// static int max_nof_steps = 0;
// static int walks_over_60 = 0;
// template <int dim>
// struct MCDriver<dim>::output MCDriver<dim>::solve(int nof_walks, const double *x_start)
// {
// 	int i, j;
// 	struct output out; //mean of the solution estimate and the nof_steps
// 	double sol_est; //temp estimate of the solution
// 	double x[dim], y[dim], d; //boundary point, ball point, sphere's radius
// 	
// 	d = calc_sphere_rad(x_start);
// 	if (d < 0.) {return out;} //outside \Omega
// 	if (d < btol) { //very close to the boundary
// 		out.msol_est = Prob.f(x_start);
// 		return out;
// 	}
// 	
// 	//out.mnof_steps = .0;
// 	out.msol_est = .0;
// 	for (i=0; i<nof_walks; i++) {
// 		for (j=0; j<dim; j++) {x[j] = x_start[j];}
// 
// 		
// 		while ((d = calc_sphere_rad(x)) > btol) {
// 			rand_update_y(x, y, d);
// 			sol_est += a(d)*Prob.q(y);
// 			
// 			rand_update_x(x, d);
// 			//out.mnof_steps++;
// // 			nof_steps++;
// 		}
// 		
// 		sol_est += Prob.f(x);
// 		
// // 		if (nof_steps > max_nof_steps) {
// // 			max_nof_steps = nof_steps;
// // 		}
// // 		if (nof_steps == 60) {
// // 			walks_over_60++;
// // 		}
// 		
// 		out.msol_est += sol_est/nof_walks;
// 	}
// 	
// 	//out.mnof_steps /= nof_walks;
// 	
// // 	std::cout<<"max_nof_steps: "<<max_nof_steps<<std::endl;
// // 	std::cout<<"walks_over_60: "<<walks_over_60<<std::endl;
// 	
// 	return out;
// }

template <int dim>
struct MCDriver<dim>::output MCDriver<dim>::solve(int nof_walks, const double *x_start)
{
	int i, j;
	struct output out; //mean of the solution estimate and the nof_steps
	double sol_est; //temp estimate of the solution
	double x[dim], y[dim], d; //boundary point, ball point, sphere's radius
	
	d = calc_sphere_rad(x_start);
	if (d < 0) {return out;} //outside \Omega
	if (d < btol) { //very close to the boundary
		out.msol_est = Prob.f(x_start);
		return out;
	}
	
#if RAND_MODE == 2
	init_quasi(5*nof_walks*30);
#endif
	
	//out.mnof_steps = .0;
	out.msol_est = .0;
	for (i=0; i<nof_walks; i++) {
		for (j=0; j<dim; j++) {x[j] = x_start[j];}
		
		sol_est = .0;
		
// 		if ((d = calc_sphere_rad(x)) > btol) {
// 			quasi_update_y(x, y, d);
// 			sol_est += a(d)*Prob.q(y);
// 			
// 			quasi_update_x(x, d);
// 			//out.mnof_steps++;
// 		}
		
		while ((d = calc_sphere_rad(x)) > btol) {
			#if RAND_MODE == 1 /*pseudo*/
				rand_update_y(x, y, d);
			#else /*quasi*/
				quasi_update_y(x, y, d);
			#endif
			sol_est += a(d)*Prob.q(y);
			
			#if RAND_MODE == 1 /*pseudo*/
				rand_update_x(x, d);
			#else /*quasi*/
				quasi_update_x(x, d);
			#endif
			//out.mnof_steps++;
		}
		
		sol_est += Prob.f(x);
		
		out.msol_est += sol_est/nof_walks;
	}
	
#if RAND_MODE == 2
	delete [] qarray;
#endif
	
	std::cout<<"diff: "<<fabs(Prob.test_u(x_start)-out.msol_est)<<std::endl;
	
	//out.mnof_steps /= nof_walks;
	return out;
}

#if INTERFACE_MODE == 2
	double D_temp[3];
#endif
template <int dim>
void* MCDriver<dim>::mcmain(void *arg)
{
	int i;
	struct mcjob<dim> mcj = *(struct mcjob<dim> *)arg;
	struct output out;
	
	for (i=0; i<mcj.nof_nodes; i++) {
#if INTERFACE_MODE == 1
		out = mcj.mc->solve(mcj.nof_walks, mcj.node_coord[i]);
		mcj.sol_est[i] = out.msol_est;
#elif INTERFACE_MODE == 2
 		Problem<dim> Prob_temp(D_temp);
		mcj.sol_est[i] = Prob_temp.test_u(mcj.node_coord[i]);
#endif
	}
	
	/* here mnof_steps could be returned */
	pthread_exit(NULL);
}

template <int dim>
void MCDriver<dim>::monte_carlo(int nof_nodes, double** node_coord, double* &sol_est)
{
#if INTERFACE_MODE == 2
	for (int tmpi=0; tmpi<dim; tmpi++) {
		D_temp[tmpi] = D[tmpi];
	}
#endif

	int i, node_offset, nof_threads;
	int nof_jpthread, ntjp; //nof jobs per thread, nof threads with plus one job
	pthread_t *pmc;
	struct mcjob<dim> *mcj;
	
	sol_est = new double[nof_nodes];
	
	nof_threads = (pp.max_nof_threads > nof_nodes) ? nof_nodes:pp.max_nof_threads;
	nof_jpthread = nof_nodes/nof_threads;
	ntjp = nof_nodes - nof_jpthread*nof_threads;
	
	std::cout<<"\n*** Monte Carlo ***"<<std::endl
	         <<"nof threads:\t"<<nof_threads<<std::endl
	         <<"nof nodes/jobs:\t"<<nof_nodes<<std::endl
	         <<"nof walks:\t"<<pp.nof_walks<<std::endl
	         <<"boundary tolerance:\t"<<pp.btol<<std::endl;
	
	pmc = new pthread_t[nof_threads];
	mcj = new mcjob<dim>[nof_threads];
	
	for (i=0; i<ntjp; i++) {
		node_offset = i*(nof_jpthread+1);
		mcj[i].mc         = this;
		mcj[i].nof_walks  = pp.nof_walks;
		mcj[i].nof_nodes  = nof_jpthread+1;
		mcj[i].node_coord = &node_coord[node_offset];
		mcj[i].sol_est    = &sol_est[node_offset];
		pthread_create(&pmc[i], NULL, mcmain, &mcj[i]);
	}
	
	for (i=ntjp; i<nof_threads; i++) {
		node_offset = ntjp*(nof_jpthread+1) + (i-ntjp)*nof_jpthread;
		mcj[i].mc         = this;
		mcj[i].nof_walks  = pp.nof_walks;
		mcj[i].nof_nodes  = nof_jpthread;
		mcj[i].node_coord = &node_coord[node_offset];
		mcj[i].sol_est    = &sol_est[node_offset];
		pthread_create(&pmc[i], NULL, mcmain, &mcj[i]);
	}
	
	for (i=0; i<nof_threads; i++) {pthread_join(pmc[i], NULL);}
// 	delete [] pmc;  delete [] mcj;
}

#endif
