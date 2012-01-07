#ifndef _PDD_H_
#define _PDD_H_

#include <iostream>
#include <ctime>
#include "dd.h"
#include "Problem.h"
#include "MCDriver.h"
#include "Inter3DDriver.h"
#include "LaplaceDriver.h"
#include "more/euclid_norm.h"

template <int dim>
class Pdd
{
private:
	int get_node_coord(const struct pdd_prm<dim> &pp, double** &nc);
public:
	void pdd(struct pdd_prm<dim> &pp);
};

//========================== 2D or 3D ==========================
#define print_mc_est(dim) \
do {std::cout<<"sol_est: "<<std::endl; \
    for (int i=0; i<nof_nodes; i++) { \
        for (int j=0; j<(dim); j++) {std::cout<<node_coord[i][j]<<" \t";} \
        std::cout<<": "<<mc_sol_est[i]<< \
                   " (test_u: "<<Prob.test_u(node_coord[i])<< \
                   ", diff:"<<(Prob.test_u(node_coord[i])-mc_sol_est[i])<<")"<< \
                   std::endl; \
    } \
} while (0)

#define print_elapsed_time() \
do {std::cout<<std::endl<<"time elapsed: "<<(time(NULL) - start_tm)<<"s"<<std::endl;} while (0)

//============================= 2D =============================
template <>
int Pdd<2>::get_node_coord(const struct pdd_prm<2> &pp, double** &nc)
{
	int i, j, nof_nodes;
	double *D;
	
	D = pp.D;
	int  Xlines = pp.dp.Xlines,  nodes_X = pp.dp.nodes_X;
	int  Ylines = pp.dp.Ylines,  nodes_Y = pp.dp.nodes_Y;
	double *Ycoord = pp.dp.Ycoord,  *Xcoord = pp.dp.Xcoord;
	
	nof_nodes = Ylines*nodes_Y + Xlines*nodes_X;
	nc = new double*[nof_nodes];
	for (i=0; i<nof_nodes; i++) {nc[i] = new double[2/*dim*/];}
	
	assert(Ylines == 0 || nodes_Y >= 2);
	assert(Xlines == 0 || nodes_X >= 2);
	
	int index = 0;
	for (i=0; i<Ylines; i++) { //1st: X coordinates
		for (j=0; j<nodes_Y; j++) {
			//int index = i*nodes_Y + j;
			nc[index][Xd] = Ycoord[i];
			nc[index][Yd] = j*(D[Yd]/(nodes_Y - 1));
			index++;
		}
	}
	
	//int allnodes_Y = Ylines*nodes_Y;
	for (i=0; i<Xlines; i++) { //2nd: Y coordinates
		for (j=0; j<nodes_X; j++) {
			//int index = allnodes_Y + j*nodes_X + j;
			nc[index][Xd] = j*(D[Xd]/(nodes_X - 1));
			nc[index][Yd] = Xcoord[i];
			index++;
		}
	}
	
	return nof_nodes;
}

//============================= 3D =============================
template <>
int Pdd<3>::get_node_coord(const struct pdd_prm<3> &pp, double** &nc)
{
	int i, j, k, nof_nodes;
	double *D;
	
	D = pp.D;
	int XYplanes = pp.dp.XYplanes[0];  int Xnodes_XY = pp.dp.XYplanes[1];  int Ynodes_XY = pp.dp.XYplanes[2];
	int YZplanes = pp.dp.YZplanes[0];  int Ynodes_YZ = pp.dp.YZplanes[1];  int Znodes_YZ = pp.dp.YZplanes[2];
	int XZplanes = pp.dp.XZplanes[0];  int Xnodes_XZ = pp.dp.XZplanes[1];  int Znodes_XZ = pp.dp.XZplanes[2];
	double *XYcoord = pp.dp.XYcoord,  *YZcoord = pp.dp.YZcoord,  *XZcoord = pp.dp.XZcoord;
	
	nof_nodes = XYplanes*Xnodes_XY*Ynodes_XY + YZplanes*Ynodes_YZ*Znodes_YZ + XZplanes*Xnodes_XZ*Znodes_XZ;
	nc = new double*[nof_nodes];
	for (i=0; i<nof_nodes; i++) {nc[i] = new double[3/*dim*/];}
	
	assert(XYplanes == 0 || Xnodes_XY >= 2);
	assert(XYplanes == 0 || Ynodes_XY >= 2);
	assert(YZplanes == 0 || Ynodes_YZ >= 2);
	assert(YZplanes == 0 || Znodes_YZ >= 2);
	assert(XZplanes == 0 || Xnodes_XZ >= 2);
	assert(XZplanes == 0 || Znodes_XZ >= 2);
	
	for (i=0; i<YZplanes; i++) { //for every YZ plane (1st: X coordinates)
		for (j=0; j<Znodes_YZ; j++) {
			for (k=0; k<Ynodes_YZ; k++) {
				int index = i*Znodes_YZ*Ynodes_YZ + j*Ynodes_YZ + k;
				nc[index][Xd] = YZcoord[i];
				nc[index][Yd] = j*(D[Yd]/(Ynodes_YZ - 1));
				nc[index][Zd] = k*(D[Zd]/(Znodes_YZ - 1));
			}
		}
	}
	
	int allYZnodes = YZplanes*Ynodes_YZ*Znodes_YZ;
	for (i=0; i<XZplanes; i++) { //for every XZ plane (2nd: Y coordinates)
		for (j=0; j<Znodes_XZ; j++) {
			for (k=0; k<Xnodes_XZ; k++) {
				int index = allYZnodes + i*Znodes_XZ*Xnodes_XZ + j*Xnodes_XZ + k;
				nc[index][Xd] = j*(D[Yd]/(Xnodes_XZ - 1));
				nc[index][Yd] = XZcoord[i];
				nc[index][Zd] = k*(D[Zd]/(Znodes_XZ - 1));
			}
		}
	}
	
	int allYZnXZnodes = allYZnodes + XZplanes*Xnodes_XZ*Znodes_XZ;
	for (i=0; i<XYplanes; i++) { //for every XY plane (3rd: Z coordinates)
		for (j=0; j<Ynodes_XY; j++) {
			for (k=0; k<Xnodes_XY; k++) {
				int index = allYZnXZnodes + i*Ynodes_XY*Xnodes_XY + j*Xnodes_XY + k;
				nc[index][Xd] = j*(D[Xd]/(Xnodes_XY - 1));
				nc[index][Yd] = k*(D[Yd]/(Ynodes_XY - 1));
				nc[index][Zd] = XYcoord[i];
			}
		}
	}
	
	return nof_nodes;
}

//========================== 2D or 3D ==========================
template <int dim>
void Pdd<dim>::pdd(struct pdd_prm<dim> &pp)
{
	time_t start_tm = time(NULL);
	
	// ==NODE'S COORDINATES==
	int nof_nodes;
	double **node_coord; //nodes' coordinates
	nof_nodes = get_node_coord(pp, node_coord); /*! IMPORTANT: get_node_coord() returns the boundary nodes also */
	
	// ==MONTE CARLO==
	double *mc_sol_est; //monte carlo solution estimates
	Problem<dim> Prob(pp.D);
	MCDriver<dim> mc(Prob, pp);
	if (nof_nodes != 0) {
		mc.monte_carlo(nof_nodes, node_coord, mc_sol_est);
		print_elapsed_time();
#define PDD_TEST 0
#if PDD_TEST == 1
		print_mc_est(dim);
		//calc_norm<dim>(node_coord, mc_sol_est, pp.D, nof_nodes);
#endif
	}
	
	/* The use of templates here may seem puzzling. But I guess noone's hurt. */
	if (dim == 2) {
		// ==INTERPOLATION==
		//2D interpolation will be set to use directly through the finite element solver (LaplaceSolve.h)
		
		// ==DETERMINISTIC==
		LaplaceDriver<dim> ld(Prob, pp);
		ld.laplace_driver(node_coord, mc_sol_est);
		print_elapsed_time();
	}
	else { // dim == 3
		// ==INTERPOLATION==
		UCBspl::SplineSurface *dec_pl; //decomposition planes
		Inter3DDriver inter3D(Prob, pp);
		if (nof_nodes != 0) {
			inter3D.interpolate(node_coord, mc_sol_est, dec_pl);
			print_elapsed_time();
		}
		
		// ==DETERMINISTIC==
		LaplaceDriver<dim> ld(Prob, pp);
		ld.laplace_driver(dec_pl/* boundary */); //in Laplace_dd.h
		print_elapsed_time();
	}
	
	//housekeeping
// 	for (int i=0; i<=nof_nodes; i++) {delete [] node_coord[i];}
// 	delete [] node_coord;
// 	delete [] mc_sol_est;
}

#endif
