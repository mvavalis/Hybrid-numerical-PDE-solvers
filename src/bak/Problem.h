#ifndef _PROBLEM_H_
#define _PROBLEM_H_

#include <cmath>
#include "dd.h"

template <int dim>
class Problem
{
private:
public:
	Problem(const double *D);
	
	double D[dim]; //\Omega: dimension length
	
	double test_u(const double *x);
	double q(const double *y);
	double f(const double *x);
};


//========================== 2D or 3D ==========================
template <int dim>
Problem<dim>::Problem(const double *D)
{int i;  for (i=0; i<dim; i++) {this->D[i] = D[i];}}


/////////// PAPER PROBLEMS

//============================= 2D =============================
template <> inline double Problem<2>::test_u(const double *x)
{return sin(M_PI*(x[0]-1))*sinh(x[1]-1)  +  cosh(2*(x[0]-1))*cos(2*M_PI*(x[1]-1));}

template <> inline double Problem<2>::q(const double *x)
{return - (1-M_PI*M_PI) * (sin(M_PI*(x[0]-1))*sinh(x[1]-1)  +  4*cosh(2*(x[0]-1))*cos(2*M_PI*(x[1]-1)));}

template <> inline double Problem<2>::f(const double *x)
{return test_u(x); /*can get it faster*/}

//============================= 3D =============================
//./main.out 3D 1  2. 2. 2.  1 1 1 u u u  0 0  0 0  0 0  1000  1e-16  n  1000 1e-15 4
template <> inline double Problem<3>::test_u(const double *x)
{return exp(M_SQRT2*M_PI*(x[0]-1)) * sin(M_PI*(x[1]-1 + x[2]-1))  +
        (1./6.) * ((x[0]-1)*(x[0]-1)*(x[0]-1)  +  (x[1]-1)*(x[1]-1)*(x[1]-1)  +  (x[2]-1)*(x[2]-1)*(x[2]-1));}

template <> inline double Problem<3>::q(const double *x)
{return - (x[0]-1 + x[1]-1 + x[2]-1);}

template <> inline double Problem<3>::f(const double *x)
{return test_u(x); /*can get it faster*/}


/////////// OLD PROBLEMS

//============================= 2D =============================
// template <> inline double Problem<2>::test_u(const double *x)
// {return x[0]*x[1]*(D[0]-x[0])*(D[1]-x[1])*exp(x[0]+x[1]);}

/* How to get q(x) automatically:
 * -- MATLAB script:
 *    clear;  syms D0 D1 x0 x1;  g=x0*x1*(D0-x0)*(D1-x1)*exp(x0+x1);
 *    q=-(diff(g,x0,2) + diff(g,x1,2));  simple(q);
 * -- Now using regular expressions turn "x.\^2" to "\0*\0" (where \0 is the "complete match" placeholder),
 *    and then "\^2" to "". Lastly turn "x0" to "x[0]", "x1" to "x[1]", "D0" to "D[0]" and "D1" to "D[1]".
 */
// template <> inline double Problem<2>::q(const double *x)
// {return -2 * exp(x[0]+x[1]) * (-x[1]*D[1] + x[1]*x[1] + x[1]*D[0]*D[1] - x[1]*x[1]*D[0] - 2*x[1]*x[0]*D[1] + 2*x[1]*x[1]*x[0] + x[0]*x[1]*D[0]*D[1] - x[0]*x[1]*x[1]*D[0] - x[0]*x[0]*x[1]*D[1] + x[0]*x[0]*x[1]*x[1] - x[0]*D[0] + x[0]*x[0] + x[0]*D[0]*D[1] - 2*x[0]*D[0]*x[1] - x[0]*x[0]*D[1] + 2*x[0]*x[0]*x[1]);}

// template <> inline double Problem<2>::test_u(const double *x)
// {return x[0]*x[1]*(D[0]-x[0])*(D[1]-x[1]);}
// 
// template <> inline double Problem<2>::q(const double *x)
// {return 2 * (  x[0]*(D[0]-x[0]) + x[1]*(D[1]-x[1])  );}
// 
// template <> inline double Problem<2>::f(const double *x)
// {return .0;}

//============================= 3D =============================
// template <> inline double Problem<3>::test_u(const double *x)
// {return x[0]*x[1]*x[2]*(D[0]-x[0])*(D[1]-x[1])*(D[2]-x[2])*exp(x[0]+x[1]+x[2]);}
// 
// template <> inline double Problem<3>::q(const double *x)
// {
// return 2*x[1]*x[2]*(D[1]-x[1])*(D[2]-x[2])*exp(x[0]+x[1]+x[2])-2*x[1]*x[2]*(D[0]-x[0])*(D[1]-x[1])*(D[2]-x[2])*exp(x[0]+x[1]+x[2])+2*x[0]*x[1]*x[2]*(D[1]-x[1])*(D[2]-x[2])*exp(x[0]+x[1]+x[2])-3*x[0]*x[1]*x[2]*(D[0]-x[0])*(D[1]-x[1])*(D[2]-x[2])*exp(x[0]+x[1]+x[2])+2*x[0]*x[2]*(D[0]-x[0])*(D[2]-x[2])*exp(x[0]+x[1]+x[2])-2*x[0]*x[2]*(D[0]-x[0])*(D[1]-x[1])*(D[2]-x[2])*exp(x[0]+x[1]+x[2])+2*x[0]*x[1]*x[2]*(D[0]-x[0])*(D[2]-x[2])*exp(x[0]+x[1]+x[2])+2*x[0]*x[1]*(D[0]-x[0])*(D[1]-x[1])*exp(x[0]+x[1]+x[2])-2*x[0]*x[1]*(D[0]-x[0])*(D[1]-x[1])*(D[2]-x[2])*exp(x[0]+x[1]+x[2])+2*x[0]*x[1]*x[2]*(D[0]-x[0])*(D[1]-x[1])*exp(x[0]+x[1]+x[2]);
// }
// 
// template <> inline double Problem<3>::f(const double *x)
// {return .0;}

#endif
