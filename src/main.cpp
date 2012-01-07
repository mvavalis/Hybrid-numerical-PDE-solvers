/**
 * More input params:
 * -Should process_solution(int) in LaplaceSolve.run() be called at the end (to get norms and stuff)?
 *  If yes a valid test_u() function in Problem.h is needed, else we don't care about test_u().
 */

#include <iostream>
#include <string.h>
#include <errno.h>
using namespace std;
#include "dd.h"
#include "Pdd.h"

void print_usage(char *argv[])
{
	cerr<<"Usage: "<<argv[0]<<" <2D | 3D> <maxt> <dlen.> <subd.> <dect.> <decc.> <nppp.> <wlks> <btol> <sels> <selsopt> <solp>"<<endl<<endl;
	
	cerr<<"maxt: max number of threads"<<endl\
		<<"dlen.: length of the domain along each dimension (2D: 2 args, 3D: 3 args)"<<endl\
		<<"subd.: number of subdomains along each dimension (2D: 2 args, 3D: 3 args)"<<endl\
		<<"dect.: decomposition type (2D: 2 args, 3D: 3 args)"<<endl\
		<<"       decomposition types are: uniform (arg: u) and non-uniform (arg: n)"<<endl\
		<<"decc.: the coordinates of the interfaces corresponding to each"<<endl\
		<<"       dimension along which the domain is decomposed non-uniformly"<<endl\
		<<"       (note: the number of interfaces along a dimension is one"<<endl\
		<<"       less than the subdomains along a dimension)"<<endl\
		<<"nppp.: number of nodes on an interface along a dimension"<<endl\
		<<"       2D: 2 args corresponding to the decomposing lines along dimension Y and X respectively"<<endl\
		<<"       (i.e. parallel to Y and parallel to X respectively)"<<endl\
		<<"       3D: 6 args, with each pair corresponding to the decomposing planes YZ, XZ and XY respectively"<<endl\
		<<"       (i.e. orthogonal to X, orthogonal to Y and orthogonal to Z respectively)"<<endl\
		<<"wlks: number of walks"<<endl\
		<<"btol: boundary tolerance"<<endl\
		<<"sels: y | n (selective solution (i.e. solve for certain subdomains only))"<<endl\
		<<"selsopt: (give this arg only if you use selective solution)"<<endl\
		<<"\tthe subdomains for which a solution is wanted; give the required input as follows:"<<endl\
		<<"\tfirstly, give the number of the subdomains, then give a list of pairs (if 2D) or"<<endl\
		<<"\ttriplets (if 3D) corresponding to the high-level coordinates of the subdomains"<<endl\
		<<"\t(the part of a coordinate corresponding to a certain dimension (say D) takes"<<endl\
		<<"\tvalues in [1, number of subdomains along dimension D])"<<endl\
		<<"solp (solver parameters): the following parameters are required for each of the"<<endl\
		<<"\tsubdomains for which a solution is wanted (if selective solution is not used,"<<endl\
		<<"\tthe paremeters are to be given only once, and will be used for all subdomains):"<<endl\
		<<"\t<max_nof_iters> <tol> <refine>"<<endl\
		<<"\tmax_nof_iters: maximum number of iteration steps that the solver takes before failure"<<endl\
		<<"\ttol: the tolerance used by the solver to determine success of the iteration"<<endl\
		<<"\trefine: laplace grid refine times"<<endl;
	cerr<<endl<<"examples:"<<endl\
		<<argv[0]<<" 2D"<<" 6"<<"  1. 1."<<"  2 4"<<" n u"<<"  .3           "<<"  12 12"<<"  5000"<<" 1e-13"<<"  n"<<"  1000 1e-15 5"<<endl\
		<<argv[0]<<" 2D"<<" 6"<<"  1. 1."<<"  2 4"<<" n n"<<"  .3  .25 .5 .7"<<"  12 12"<<"  5000"<<" 1e-13"<<"  n"<<"  1000 1e-15 5"<<endl\
		<<argv[0]<<" 2D"<<" 6"<<"  1. 1."<<"  2 4"<<" n n"<<"  .3  .25 .5 .7"<<"  12 12"<<"  5000"<<" 1e-13"<<"  y"<<" 3  1 1  1 3  2 3"<<"  1000 1e-15 5  1000 1e-15 2  500 1e-10 2"<<endl\
		<<argv[0]<<" 3D"<<" 8"<<"  1. 1. 1.5"<<"  3 3 2"<<" u u u"<<"                    "<<"  5 5  5 5  7 7"<<"  10000"<<" 1e-16"<<"  n"<<"  1000 1e-15 4"<<endl\
		<<argv[0]<<" 3D"<<" 8"<<"  1. 1. 1.5"<<"  3 3 2"<<" u n u"<<"         .5 .55     "<<"  5 5  5 5  7 7"<<"  10000"<<" 1e-16"<<"  n"<<"  1000 1e-15 4"<<endl\
		<<argv[0]<<" 3D"<<" 8"<<"  1. 1. 1.5"<<"  3 3 2"<<" n n n"<<"  .3 .4  .5 .55  1.2"<<"  5 5  5 5  7 7"<<"  10000"<<" 1e-16"<<"  n"<<"  1000 1e-15 4"<<endl\
		<<argv[0]<<" 3D"<<" 8"<<"  1. 1. 1.5"<<"  3 3 2"<<" n n n"<<"  .3 .4  .5 .55  1.2"<<"  5 5  5 5  7 7"<<"  10000"<<" 1e-16"<<"  y"<<" 2  1 2 1  3 3 1"<<"  3000 1e-17 3  1000 1e-15 2"<<endl;
	exit(1);
}

#define in_mc(n) /*more check: check if there are @n more args (after @i)*/ \
do {if (i+((n)-1) > argc-1) {print_usage(argv);}} while (0)

#define read_int(to, from) /*don't use @to or @from twice in here*/ \
do { \
	in_mc(1); \
	(to) = strtol((from), (char **) NULL, 10); \
	if (errno == ERANGE) {print_usage(argv);} \
} while(0)

#define read_double(to, from) /*don't use @to or @from twice in here*/ \
do { \
	in_mc(1); \
	(to) = strtod((from), NULL); \
	if (errno == ERANGE) {print_usage(argv);} \
} while(0)

//read <subd.> <dect.> <decc.> and <nppp.>
int set_dec_prm(struct pdd_prm<2> &pp, int argc, char *argv[], int i)
{
	int j;
	bool Yunif, Xunif; //if the domain is uniformly decomposed along dimensions Y and X
	
	read_int(pp.dp.Ylines, argv[i++]);  pp.dp.Ylines--;  assert(pp.dp.Ylines>=0); //read <subd.>
	read_int(pp.dp.Xlines, argv[i++]);  pp.dp.Xlines--;  assert(pp.dp.Xlines>=0);
	cout<<"Ylines: "<<pp.dp.Ylines<<endl;
	cout<<"Xlines: "<<pp.dp.Xlines<<endl;
	in_mc(1);  Yunif = (argv[i][0] == 'u' || argv[i][0] == 'U') ? true:false;  i++; //read <dect.>
	in_mc(1);  Xunif = (argv[i][0] == 'u' || argv[i][0] == 'U') ? true:false;  i++;
	cout<<"Yunif: "<<Yunif<<endl;
	cout<<"Xunif: "<<Xunif<<endl;
	pp.dp.Ycoord = new double[pp.dp.Ylines];  pp.dp.Xcoord = new double[pp.dp.Xlines]; //read <decc.>
	for (j=0; j<pp.dp.Ylines; j++) {
		if (Yunif) {pp.dp.Ycoord[j] = (j+1)*(pp.D[Xd]/(pp.dp.Ylines + 1));}
		else {read_double(pp.dp.Ycoord[j], argv[i++]);}
	}
	for (j=0; j<pp.dp.Xlines; j++) {
		if (Xunif) {pp.dp.Xcoord[j] = (j+1)*(pp.D[Yd]/(pp.dp.Xlines + 1));}
		else {read_double(pp.dp.Xcoord[j], argv[i++]);}
	}
	read_int(pp.dp.nodes_Y, argv[i++]); //read <nppp.>
	read_int(pp.dp.nodes_X, argv[i++]);
	cout<<"nodes_Y: "<<pp.dp.nodes_Y<<endl;
	cout<<"nodes_X: "<<pp.dp.nodes_X<<endl;
		//if no interfaces, then no nodes along them
		assert((pp.dp.Ylines>0 && pp.dp.nodes_Y>0) || (pp.dp.Ylines==0 && pp.dp.nodes_Y==0));
		assert((pp.dp.Xlines>0 && pp.dp.nodes_X>0) || (pp.dp.Xlines==0 && pp.dp.nodes_X==0));
	return i;
}

//read <subd.> <dect.> <decc.> and <nppp.>
int set_dec_prm(struct pdd_prm<3> &pp, int argc, char *argv[], int i)
{
	int j;
	bool YZunif, XZunif, XYunif; //if the domain is uniformly decomposed along dimensions X, Y and Z, respectively
	
	cout<<"i: "<<i<<endl;
	read_int(pp.dp.YZplanes[0], argv[i++]);  pp.dp.YZplanes[0]--;  assert(pp.dp.YZplanes[0]>=0); //read <subd.>
	read_int(pp.dp.XZplanes[0], argv[i++]);  pp.dp.XZplanes[0]--;  assert(pp.dp.XZplanes[0]>=0);
	read_int(pp.dp.XYplanes[0], argv[i++]);  pp.dp.XYplanes[0]--;  assert(pp.dp.XYplanes[0]>=0);
	cout<<"YZplanes[0]: "<<pp.dp.YZplanes[0]<<endl;
	cout<<"XZplanes[0]: "<<pp.dp.XZplanes[0]<<endl;
	cout<<"XYplanes[0]: "<<pp.dp.XYplanes[0]<<endl;
	cout<<"i: "<<i<<endl;
	in_mc(1);  YZunif = (argv[i][0] == 'u' || argv[i][0] == 'U') ? true:false;  i++; //read <dect.>
	in_mc(1);  XZunif = (argv[i][0] == 'u' || argv[i][0] == 'U') ? true:false;  i++;
	in_mc(1);  XYunif = (argv[i][0] == 'u' || argv[i][0] == 'U') ? true:false;  i++;
	cout<<"YZunif: "<<YZunif<<endl;
	cout<<"XZunif: "<<XZunif<<endl;
	cout<<"XYunif: "<<XYunif<<endl;
	cout<<"i: "<<i<<endl;
	pp.dp.YZcoord = new double[pp.dp.YZplanes[0]]; //read <decc.>
	pp.dp.XZcoord = new double[pp.dp.XZplanes[0]];
	pp.dp.XYcoord = new double[pp.dp.XYplanes[0]];
	for (j=0; j<pp.dp.YZplanes[0]; j++) {
		if (YZunif) {pp.dp.YZcoord[j] = (j+1)*(pp.D[Xd]/(pp.dp.YZplanes[0] + 1));}
		else {read_double(pp.dp.YZcoord[j], argv[i++]);}
	}
	for (j=0; j<pp.dp.XZplanes[0]; j++) {
		if (XZunif) {pp.dp.XZcoord[j] = (j+1)*(pp.D[Yd]/(pp.dp.XZplanes[0] + 1));}
		else {read_double(pp.dp.XZcoord[j], argv[i++]);}
	}
	for (j=0; j<pp.dp.XYplanes[0]; j++) {
		if (XYunif) {pp.dp.XYcoord[j] = (j+1)*(pp.D[Zd]/(pp.dp.XYplanes[0] + 1));}
		else {read_double(pp.dp.XYcoord[j], argv[i++]);}
	}
	cout<<"i: "<<i<<endl;
	read_int(pp.dp.YZplanes[1], argv[i++]);  read_int(pp.dp.YZplanes[2], argv[i++]); //read <nppp.>
	read_int(pp.dp.XZplanes[1], argv[i++]);  read_int(pp.dp.XZplanes[2], argv[i++]);
	read_int(pp.dp.XYplanes[1], argv[i++]);  read_int(pp.dp.XYplanes[2], argv[i++]);
	cout<<"YZplanes: "<<pp.dp.YZplanes[1]<<" "<<pp.dp.YZplanes[2]<<endl;
	cout<<"XZplanes: "<<pp.dp.XZplanes[1]<<" "<<pp.dp.XZplanes[2]<<endl;
	cout<<"XYplanes: "<<pp.dp.XYplanes[1]<<" "<<pp.dp.XYplanes[2]<<endl;
	cout<<"i: "<<i<<endl;
		//if no interfaces, then no nodes along them
		assert((pp.dp.YZplanes[0]==0 && pp.dp.YZplanes[1]==0 && pp.dp.YZplanes[2]==0) ||
		       (pp.dp.YZplanes[0] >0 && pp.dp.YZplanes[1] >0 && pp.dp.YZplanes[2] >0));
		assert((pp.dp.XZplanes[0]==0 && pp.dp.XZplanes[1]==0 && pp.dp.XZplanes[2]==0) ||
		       (pp.dp.XZplanes[0] >0 && pp.dp.XZplanes[1] >0 && pp.dp.XZplanes[2] >0));
		assert((pp.dp.XYplanes[0]==0 && pp.dp.XYplanes[1]==0 && pp.dp.XYplanes[2]==0) ||
		       (pp.dp.XYplanes[0] >0 && pp.dp.XYplanes[1] >0 && pp.dp.XYplanes[2] >0));
	return i;
}

//read <sels> <selsopt> <solp>
template <int dim>
int set_solve_prm(struct pdd_prm<dim> &pp, int argc, char *argv[], int i)
{
	in_mc(1);  pp.sp.selective_sol = (argv[i][0] == 'y' || argv[i][0] == 'y') ? true:false;  i++; //read <sels>
	
	 //read <selsopt>
	if (pp.sp.selective_sol) {read_int(pp.sp.nof_subd, argv[i++]);}
	else                     {pp.sp.nof_subd = 1;}
	
	pp.sp.subp = new solve_prm::subd_solve_prm[pp.sp.nof_subd];
	if (pp.sp.selective_sol) {
		for (int j=0; j<pp.sp.nof_subd; j++) {
			pp.sp.subp[j].subd_coord = new double[dim];
			for (int k=0; k<dim; k++) {
				read_int(pp.sp.subp[j].subd_coord[k], argv[i++]);
				cout<<"subd_coord["<<k<<"]: "<<pp.sp.subp[j].subd_coord[k]<<endl;
			}
		}
	}
	
	 //read <solp>
	for (int j=0; j<pp.sp.nof_subd; j++) {
		read_int   (pp.sp.subp[j].max_nof_iters, argv[i++]);
		cout<<"max_nof_iters: "<<pp.sp.subp[j].max_nof_iters<<endl;
		read_double(pp.sp.subp[j].tol          , argv[i++]);
		cout<<"tol: "<<pp.sp.subp[j].tol<<endl;
		read_int   (pp.sp.subp[j].refine_times , argv[i++]);
		cout<<"refine_times: "<<pp.sp.subp[j].refine_times<<endl;
	}
	
	return i;
}


template <int dim>
void init(int argc, char *argv[], struct pdd_prm<dim> &pp)
{
	int i, j;
	
	//if (!((argc==12 && !strcmp(argv[i=1], "2D")) || (argc==18 && !strcmp(argv[i=1], "3D")))) {print_usage(argv);}
	
	i=2;
	read_int(pp.max_nof_threads, argv[i++]); //read <maxt>
	cout<<"max_nof_threads: "<<pp.max_nof_threads<<endl;
	pp.D = new double[dim];  for (j=0; j<dim; j++) {read_double(pp.D[j], argv[i++]);} //read <dlen.>
	cout<<"dlen.: ";  for (j=0; j<dim; j++) {cout<<pp.D[j]<<" ";}  cout<<endl;
	i = set_dec_prm(pp, argc, argv, i); //read <subd.> <dect.> <decc.> and <nppp.>
	read_int(pp.nof_walks, argv[i++]); //read <wlks>
	cout<<"nof_walks: "<<pp.nof_walks<<endl;
	read_double(pp.btol, argv[i++]); //read <btol>
	cout<<"btol: "<<pp.btol<<endl;
	i = set_solve_prm(pp, argc, argv, i); //read <sels> <selsopt> <solp>
	
	cout.precision(10);
	cout<<"##############################################"<<endl
		<<"# Solution of Poisson's equation through the #"<<endl
		<<"# use of Probabilistic Domain Decomposition: #"<<endl
		<<"# - Monte Carlo                              #"<<endl
		<<"# - Splines                                  #"<<endl
		<<"# - Finite Elements                          #"<<endl
		<<"##############################################"<<endl;
}

int main(int argc, char *argv[])
{
	int dim;
	
	if (argc<2 || (strcmp(argv[1], "2D") && strcmp(argv[1], "3D"))) {print_usage(argv);}
	else {dim = !strcmp(argv[1], "2D") ? 2:3;}
	
	if (dim == 2) {
		struct pdd_prm<2> pp;
		init(argc, argv, pp);
		Pdd<2> pdd;  pdd.pdd(pp);
	}
	else {
		struct pdd_prm<3> pp;
		init(argc, argv, pp);
		Pdd<3> pdd;  pdd.pdd(pp);
	}
	
	return 0;
}
