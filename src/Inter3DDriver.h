#ifndef _INTER_H_
#define _INTER_H_

#include <iostream>
#include <MBA.h>
#include <UCButils.h>
#include "dd.h"
#include "Pdd.h"
#include "Problem.h"

//============================= 3D =============================
class Inter3DDriver
{
private:
	Problem<3> &Prob;
	const struct pdd_prm<3> &pp;
	
	static void* intermain(void *arg);
public:
	Inter3DDriver(Problem<2> &Prob, const struct pdd_prm<2> &pp); //this is a dummy function
	Inter3DDriver(Problem<3> &Prob, const struct pdd_prm<3> &pp);
	void interpolate(double **node_coord, double *mc_sol_est, UCBspl::SplineSurface* &dec_pl);
};

struct interdata {
	boost::shared_ptr<dVec> *x_arr;
	boost::shared_ptr<dVec> *y_arr;
	boost::shared_ptr<dVec> *z_arr;
	UCBspl::SplineSurface *dec_pl;
};

struct interjob {
	struct interdata *intdata;
	int nof_planes;
	int plane_offset;
};

Inter3DDriver::Inter3DDriver(Problem<3> &Prob, const struct pdd_prm<3> &pp) : Prob(Prob), pp(pp)
{}

void* Inter3DDriver::intermain(void *arg)
{
	struct interjob intj = *(struct interjob *)arg;
	
	for (int i=intj.plane_offset; i<intj.plane_offset+intj.nof_planes; i++) {
//		/** TEMP (use real solution instead of estimated) */
// 		double ppD[3] = {2., 2., 2.};
// 		Problem<3> Probb(ppD);
// 		for (unsigned int j=0; j<intj.intdata->x_arr[i]->size(); j++) {
// 			dVec d1 = *intj.intdata->x_arr[i];
// 			dVec d2 = *intj.intdata->y_arr[i];
// 			double x[3] = {1., d1[j], d2[j]};
// 			(*intj.intdata->z_arr[i])[j] = Probb.test_u(x);
// 		}
//		/** ~TEMP */
		
		MBA mba(intj.intdata->x_arr[i], intj.intdata->y_arr[i], intj.intdata->z_arr[i]);
		
		/* set the first two params of MBAalg() according to the dimensions of the plane */
		double usize = *std::max_element(intj.intdata->x_arr[i]->begin(), intj.intdata->x_arr[i]->end());
		double vsize = *std::max_element(intj.intdata->y_arr[i]->begin(), intj.intdata->y_arr[i]->end());
		int m0, n0;
		if (usize > vsize) {m0 = (int)((usize/vsize)+.5);  n0 = 1;}
		else               {m0 = 1;  n0 = (int)((vsize/usize)+.5);}
		mba.MBAalg(m0,n0,10);
		std::cout<<"m0: "<<m0<<std::endl;
		std::cout<<"n0: "<<n0<<std::endl;
		//mba.MBAalg(1,1,7);
		
		intj.intdata->dec_pl[i] = mba.getSplineSurface();
	std::cout<<"Inter3DDriver(): umax: "<<intj.intdata->dec_pl[i].umax()<<std::endl;
	std::cout<<"         umin: "<<intj.intdata->dec_pl[i].umin()<<std::endl;
	std::cout<<"         vmax: "<<intj.intdata->dec_pl[i].vmax()<<std::endl;
	std::cout<<"         vmin: "<<intj.intdata->dec_pl[i].vmin()<<std::endl;
	}
	
	pthread_exit(NULL);
}

void Inter3DDriver::interpolate(double **node_coord, double *mc_sol_est, UCBspl::SplineSurface* &dec_pl)
{
	int i, j, l, nof_threads, nof_planes;
	int nof_jpthread, ntjp; //nof jobs per thread, nof threads with plus one job
	pthread_t *pmc;
	
	int XYplanes = pp.dp.XYplanes[0];  int Xnodes_XY = pp.dp.XYplanes[1];  int Ynodes_XY = pp.dp.XYplanes[2];
	int YZplanes = pp.dp.YZplanes[0];  int Ynodes_YZ = pp.dp.YZplanes[1];  int Znodes_YZ = pp.dp.YZplanes[2];
	int XZplanes = pp.dp.XZplanes[0];  int Xnodes_XZ = pp.dp.XZplanes[1];  int Znodes_XZ = pp.dp.XZplanes[2];
	
	dec_pl = new UCBspl::SplineSurface[YZplanes + XZplanes + XYplanes];
	
	struct interjob *intj;
	double tmp_point[3/*dim*/-1];
	
	nof_planes = YZplanes + XZplanes + XYplanes;
	
	typedef std::vector<double> dVec;
	boost::shared_ptr<dVec> x_arr[nof_planes];
	boost::shared_ptr<dVec> y_arr[nof_planes];
	boost::shared_ptr<dVec> z_arr[nof_planes];
	for (i=0; i<nof_planes; i++) {
		x_arr[i] = boost::shared_ptr<dVec>(new dVec);
		y_arr[i] = boost::shared_ptr<dVec>(new dVec);
		z_arr[i] = boost::shared_ptr<dVec>(new dVec);
	}
	
#define print_point(pos) \
do {std::cout<<(*x_arr[(pos)])[x_arr[(pos)]->size()-1] \
		<<" \t"<<(*y_arr[(pos)])[y_arr[(pos)]->size()-1] \
		<<" \t:"<<(*z_arr[(pos)])[z_arr[(pos)]->size()-1]<<std::endl; \
} while (0)
	
#define push_point(pos, z_val) \
do {(*x_arr[(pos)]).push_back(tmp_point[0]); \
    (*y_arr[(pos)]).push_back(tmp_point[1]); \
    (*z_arr[(pos)]).push_back((z_val)); \
				print_point((pos)); \
} while (0)
	
	std::cout<<"*** Interpolation ***"<<std::endl;
	
	std::cout<<"YZplanes:"<<YZplanes<<", Ynodes_YZ:"<<Ynodes_YZ<<", Znodes_YZ:"<<Znodes_YZ<<std::endl;
	for (l=0; l<YZplanes; l++) { //for every YZplane
		for (i=0; i<Ynodes_YZ; i++) {
			for (j=0; j<Znodes_YZ; j++) {
				int index = i*Znodes_YZ + j;
				tmp_point[0] = node_coord[index][Yd];
				tmp_point[1] = node_coord[index][Zd];
				push_point(l,  mc_sol_est[index]);
			}
		}
		//dims[l]->push_back(pp.D[Yd]);  dims[l]->push_back(pp.D[Zd]);
	}
	
	std::cout<<"XZplanes:"<<XZplanes<<", Xnodes_XZ:"<<Xnodes_XZ<<", Znodes_XZ:"<<Znodes_XZ<<std::endl;
	for (l=0; l<XZplanes; l++) { //for every XZplane
		int plane_vec_offset = YZplanes;
		int node_offset = YZplanes*Ynodes_YZ*Znodes_YZ;
		
		for (i=0; i<Xnodes_XZ; i++) {
			for (j=0; j<Znodes_XZ; j++) {
				int index = node_offset + i*Znodes_XZ + j;
				tmp_point[0] = node_coord[index][Xd];
				tmp_point[1] = node_coord[index][Zd];
				push_point(plane_vec_offset + l, mc_sol_est[index]);
			}
		}
		//dims[l]->push_back(pp.D[Xd]);  dims[l]->push_back(pp.D[Zd]);
	}
	
	std::cout<<"XYplanes:"<<XYplanes<<", Xnodes_XY:"<<Xnodes_XY<<", Ynodes_XY:"<<Ynodes_XY<<std::endl;
	for (l=0; l<XYplanes; l++) { //for every XYplane
		int plane_vec_offset = XZplanes + YZplanes;
		int node_offset = XZplanes*Xnodes_XZ*Znodes_XZ + YZplanes*Ynodes_YZ*Znodes_YZ;
		
		for (i=0; i<Xnodes_XY; i++) {
			for (j=0; j<Ynodes_XY; j++) {
				int index = node_offset + i*Ynodes_XY + j;
				tmp_point[0] = node_coord[index][Xd];
				tmp_point[1] = node_coord[index][Yd];
				push_point(plane_vec_offset + l, mc_sol_est[index]);
			}
		}
		//dims[l]->push_back(pp.D[Xd]);  dims[l]->push_back(pp.D[Yd]);
	}
	
	nof_threads = (pp.max_nof_threads > nof_planes) ? nof_planes:pp.max_nof_threads;
	nof_jpthread = nof_planes/nof_threads;
	ntjp = nof_planes - nof_jpthread*nof_threads;
	
	pmc = new pthread_t[nof_threads];
	intj = new interjob[nof_threads];
	struct interdata intdata = {x_arr, y_arr, z_arr, /*dims,*/ dec_pl};
	
	for (i=0; i<ntjp; i++) {
		intj[i].plane_offset = i*(nof_jpthread+1);
		intj[i].nof_planes   = nof_jpthread+1;
		intj[i].intdata      = &intdata;
		pthread_create(&pmc[i], NULL, intermain, &intj[i]);
	}
	for (i=ntjp; i<nof_threads; i++) {
		intj[i].plane_offset = ntjp*(nof_jpthread+1) + (i-ntjp)*nof_jpthread;
		intj[i].nof_planes   = nof_jpthread;
		intj[i].intdata      = &intdata;
		pthread_create(&pmc[i], NULL, intermain, &intj[i]);
	}
	
	for (i=0; i<nof_threads; i++) {pthread_join(pmc[i], NULL);}
	//delete [] pmc;  delete [] intj;
}

//====================== dummy functions =======================
Inter3DDriver::Inter3DDriver(Problem<2> &Prob, const struct pdd_prm<2> &pp) : Prob(*(Problem<3> *)NULL), pp(*(pdd_prm<3> *)NULL)
{}

#endif
