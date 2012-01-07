#ifndef _DD_H_
#define _DD_H_

//=============Monte Carlo Parameters================
#define RAND_MODE 1 /* 1: pseudo-rand, 2: quasi-rand */
//===================================================

//=================interface mode====================
/**
 * 1: use interpolated monte carlo values
 * 2: use interpolated real values
 * 3: use real values
 */
#define INTERFACE_MODE 1
#define TEST_INTERFACE 1 /* works for 2D only */
//===================================================

//=============Interpolation Parameters==============
#define INTERP2D_LIB 2 /* 1: Burkardt, 2: GSL */
//===================================================

//================some handy enums===================
enum {Xd=0, Yd, Zd}; //dimensions
enum {XYp=0, YZp, XZp}; //planes
//===================================================

//=============Decomposition Parameters==============
template <int dim>
struct dec_prm;

template <>
struct dec_prm<2> {
	int Ylines, Xlines; //nof lines parallel to y (horizontal), x (vertical) axis respectively
	int nodes_Y, nodes_X; //nof nodes along Y, X dimensions respectively (corresponding to Ylines, Xlines respectively)
	double *Ycoord, *Xcoord; //the coordinates of the interfaces corresponding to Ylines and Xlines, respectively
	                         //(arrays of size @Ylines and @Xlines respectively)
};

template <>
struct dec_prm<3> {
	//[0]: nof planes. [1], [2]: nof nodes along 1st and 2nd {dimension} on the {plane} respectively
	int YZplanes[3], XZplanes[3], XYplanes[3];
	double *YZcoord, *XZcoord, *XYcoord;
};
//===================================================

//=============Laplace Solver Parameters=============
struct solve_prm {
	struct subd_solve_prm {
		double *subd_coord; //the coordinates of the subdomain (used only if selective_sol==true)
		
		unsigned int max_nof_iters; //maximum number of iteration steps before failure
		double tol; //tolerance to determine success of the iteration
		unsigned int refine_times; //the nof times to refine the triangulation's cells
	};
	
	bool selective_sol; //solve for certain subdomains only
	
	int nof_subd; //nof subdomains (if (selective_sol==false) {nof_subd = 1;})
	struct subd_solve_prm *subp; //subdomains' parameters (of size @nof_subd)
};
//===================================================

//=================Pdd Parameters====================
template <int dim>
struct pdd_prm {
	//general
	int max_nof_threads;
	double *D; //\Omega: dimension length
	
	//monte carlo
	int nof_walks;
	double btol;
	
	//interpolation
	struct dec_prm<dim> dp;
	
	//laplace
	struct solve_prm sp;
};
//===================================================

#endif
