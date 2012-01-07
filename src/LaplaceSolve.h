/**
 * Based on the tutorials in the documentation of deal.II
 *   library's version 6.1.0.
 * The original can be found on deal.II's homepage
 *   (currently http://www.dealii.org)
 */

/**
 * Author: Wolfgang Bangerth, University of Heidelberg, 1999
 *
 *    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005,
 *                  2006, 2007, 2008 by the deal.II authors
 *
 *    This file is subject to QPL and may not be  distributed
 *    without copyright and license information. Please refer
 *    to the file deal.II/doc/license.html for the  text  and
 *    further information on this license.
 */

#ifndef _LAPLACESOLVE_H_
#define _LAPLACESOLVE_H_

#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_q.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>

#include <numerics/data_out.h>
#include <fstream>
#include <iostream>

#include <base/logstream.h>

#include <string>
#include <sys/stat.h>

#include <base/convergence_table.h>

#include "Problem.h"
#if   INTERP2D_LIB == 1
	#include "Burkardt/spline/spline.H"
#elif INTERP2D_LIB == 2
	#include <gsl/gsl_errno.h>
	#include <gsl/gsl_spline.h>
#endif
#include <MBA.h>
#include <UCButils.h>

using namespace dealii;

//=============Decomposition Hyperplane==============
template <int dim> struct dec_hpl;
template <> struct dec_hpl<2> {
	int nof_nodes;
#if   INTERP2D_LIB == 1
	double *spl_node_coord;
	double *mc_sol_est;
#elif INTERP2D_LIB == 2
	gsl_interp_accel *acc;
	gsl_spline *spline;
#endif
};
template <> struct dec_hpl<3> {
	UCBspl::SplineSurface *dec_pl;
};
//===================================================

//============Laplace Problem Parameters=============
template <int dim>
struct lapprm {
	//==============data==============
	Problem<dim> &Prob;
	struct dec_hpl<dim> *hpl;
	//================================
	
	//=======data handling info=======
	std::vector<Point<dim> > &grid_hrect_p1, &grid_hrect_p2; //hyper rectangle
	
	//if the left or right (init_bound[job_id][i][0] or init_bound[job_id][i][1] respectively)
	//bound of dimension i corresponds to the initial boundary of the problem
	//(if not, it corresponds to one of the decomposition hyperplanes).
	bool ***init_bound; //[job_id][dim][2]
	
	//if (!init_bound[job_id][i][j]) {
	//  2D: dim1: for (k=0; k<nof_nodes; k++) {hpl[hp_spline[job_id][i][j]].spl_node_coord[k]}, and
	//      dim2: for (k=0; k<nof_nodes; k++) {hpl[hp_spline[job_id][i][j]].mc_sol_est[k]}
	//      will be interpolated in order to acquire the boundary
	//  3D: hpl[hp_spline[i][j]].dec_pl is the boundary
	//}
	int ***hp_spline; //[job_id][dim][2], hyperplane spline
	//================================
	
	//=============params=============
	const struct solve_prm *sp;
	//================================
};
//===================================================

//=====================Solution======================
template <int dim>
class Solution : public Function<dim>
{
	private:
		Problem<dim> &Prob;
	public:
		Solution(Problem<dim> &Prob) : Function<dim>(), Prob(Prob) {}
		virtual double value          (const Point<dim> &p, const unsigned int  component = 0) const
		{
			if (dim==2) {double x[2] = {p[0], p[1]};        return Prob.test_u(x);}
			if (dim==3) {double x[3] = {p[0], p[1], p[2]};  return Prob.test_u(x);}
		}
// 		virtual Tensor<1,dim> gradient(const Point<dim> &p, const unsigned int  component = 0) const
// 		{
// 			double x[2] = {p[0], p[1]};
// 			
// 			Tensor<1,dim> return_value;
// 			if (dim==2) {
// 				return_value[0] = M_PI*cos(M_PI*(x[0]-1))*sinh(x[1]-1)  +  2*sinh(2*(x[0]-1))*cos(2*M_PI*(x[1]-1));
// 				return_value[1] = sin(M_PI*(x[0]-1))*cosh(x[1]-1)  -  2*M_PI*cosh(2*(x[0]-1))*sin(2*M_PI*(x[1]-1));
// 			}
// 			return return_value;
// 		}
};
//===================================================

//==================Right Hand Side==================
template <int dim> class LaplaceSolve;

double temp;
void *rhs;
template <int dim>
class RightHandSide : public Function<dim>
{
private:
	Problem<dim> &Prob;
public:
	RightHandSide(Problem<dim> &Prob) : Function<dim>(), Prob(Prob) {}
	
	//2D
	virtual double value(double d1, double d2,            const unsigned int component = 0) const
	{double x[2] = {d1, d2};      return Prob.q(x);}
	//3D
	virtual double value(double d1, double d2, double d3, const unsigned int component = 0) const
	{double x[3] = {d1, d2, d3};  return Prob.q(x);}
};
//===================================================

//==================Boundary Values==================
template <int dim>
class BoundaryValues : public Function<dim>
{
private:
	int job_id;
	const struct lapprm<dim> &lapp;
	
	//new boundary value (created through monte carlo estimation and interpolation)
	virtual double non_init_bound_value(const Point<dim> &p, struct dec_hpl<dim> *hpl, int min_dim) const;
public:
	BoundaryValues(int job_id, struct lapprm<dim> &lapp) : Function<dim>(), lapp(lapp) {this->job_id = job_id;}
	virtual double value(const Point<dim> &p, const unsigned int component = 0) const;
};

//==== 2D ====
template <>
double BoundaryValues<2>::non_init_bound_value(const Point<2> &p, struct dec_hpl<2> *hpl, int min_dim) const
{
	assert(hpl->nof_nodes);
#if   INTERP2D_LIB == 1
	return spline_b_val(hpl->nof_nodes, hpl->spl_node_coord, hpl->mc_sol_est, (min_dim==Xd) ? p[1]:p[0]);
#elif INTERP2D_LIB == 2
	return gsl_spline_eval(hpl->spline, (min_dim==Xd) ? p[1]:p[0], hpl->acc);
#endif
}

//==== 3D ====
template <>
double BoundaryValues<3>::non_init_bound_value(const Point<3> &p, struct dec_hpl<3> *hpl, int min_dim) const
{return (*hpl->dec_pl).f((min_dim==Xd)?p[1]:p[0], (min_dim==Zd)?p[1]:p[2]);}

//== 2D or 3D ==
template <int dim>
double BoundaryValues<dim>::value(const Point<dim> &p, const unsigned int) const
{
	int i;
	bool left; //left (or right?) bound
	int min_dim; // Xd, Yd, Zd
	double min, tmp;
	
// 	std::cout<<"BoundaryValues::value(). p(i): ";
// 	for (i=0; i<dim; i++) {
// 		std::cout<<p(i);
// 		if (i != dim-1) {std::cout<<"\t";}
// 	}
	
	//among the hyperplanes, of which the hyper rectangle of the problem consists, find the one closest to the point @p
	min = p[Xd] - lapp.grid_hrect_p1[job_id][Xd];  min_dim = Xd;  left = true;
	for (i=1; i<dim; i++) {if ((tmp = p[i] - lapp.grid_hrect_p1[job_id][i]) < min) {min = tmp;  min_dim = i;               }}
	for (i=0; i<dim; i++) {if ((tmp = lapp.grid_hrect_p2[job_id][i] - p[i]) < min) {min = tmp;  min_dim = i;  left = false;}}
	
	if (lapp.init_bound[job_id][min_dim][left ? 0:1]) { //initial boundary
		double *x = new double[dim];  for (int i=0; i<dim; i++) {x[i]=p[i];}
		double bv = lapp.Prob.f(x);
		/*delete [] x;*/  return bv;
	}
	else { //new boundary (estimated)
#if   INTERFACE_MODE != 3
		struct dec_hpl<dim> *hpl = &lapp.hpl[lapp.hp_spline[job_id][min_dim][left ? 0:1]];
		return non_init_bound_value(p, hpl, min_dim);
#elif INTERFACE_MODE == 3
 		double *x = new double[dim];  for (int i=0; i<dim; i++) {x[i]=p[i];}
 		double bv = lapp.Prob.test_u(x);
 		return bv;
#endif
// 		double bv_est = non_init_bound_value(p, hpl, min_dim);
// 		
// 		double *x = new double[dim];  for (int i=0; i<dim; i++) {x[i]=p[i];}
// 		double bv = lapp.Prob.test_u(x);
// // 		std::cout<<x[0]<<"\t"<<x[1]<<"\t"<<x[2]<<"\t: bv: "<<bv<<",\tbv_est: "<<bv_est<<",\tdiff: "<<(bv-bv_est)<<"\n";
// 		
// 		double diff[9];
// 		diff[0] = fabs(bv - non_init_bound_value(p, &lapp.hpl[0], 0));
// 		diff[1] = fabs(bv - non_init_bound_value(p, &lapp.hpl[0], 1));
// 		diff[2] = fabs(bv - non_init_bound_value(p, &lapp.hpl[0], 2));
// 		diff[3] = fabs(bv - non_init_bound_value(p, &lapp.hpl[1], 0));
// 		diff[4] = fabs(bv - non_init_bound_value(p, &lapp.hpl[1], 1));
// 		diff[5] = fabs(bv - non_init_bound_value(p, &lapp.hpl[1], 2));
// 		diff[6] = fabs(bv - non_init_bound_value(p, &lapp.hpl[2], 0));
// 		diff[7] = fabs(bv - non_init_bound_value(p, &lapp.hpl[2], 1));
// 		diff[8] = fabs(bv - non_init_bound_value(p, &lapp.hpl[2], 2));
// 		double diffmin = diff[0];
// 		for (int i=1; i<9; i++) {
// 			if (diff[i] < diffmin) {diffmin = diff[i];}
// 		}
// 		if (diffmin < fabs(bv-bv_est)) {
// // 			std::cout<<x[0]<<"\t"<<x[1]<<"\t"<<x[2]<<"\t: bv: "<<bv<<",\tbv_est: "<<bv_est<<",\tdiff: "<<(bv-bv_est)<<"\n";
// // 			std::cout<<"\t\t"<<"diffmin: "<<diffmin<<"\n";
// 			std::cout<<"\t\t"<<"relative fault: "<<fabs(bv-bv_est)/diffmin<<"\n";
// 		}
// 		
// 		return bv_est;
	}
}
//===================================================

//==================Laplace Problem==================
template <int dim>
class LaplaceSolve
{
public:
	LaplaceSolve(int basis_poly_deg, int job_id, struct lapprm<dim> &lapp)
		: fe(basis_poly_deg), dof_handler(triangulation), lapp(lapp) {
		this->job_id = job_id;
		this->grid_hrect_p1 = &lapp.grid_hrect_p1[job_id];
		this->grid_hrect_p2 = &lapp.grid_hrect_p2[job_id];
		this->basis_poly_deg = basis_poly_deg;
	}
	void run();
	Vector<double> get_solution() {return solution;}
	void process_solution(const unsigned int cycle);
	
private:
	void setup_system();
	void assemble_system();
	void solve();
	void output_results() const;
	
	Triangulation<dim>   triangulation;
	FE_Q<dim>            fe;
	DoFHandler<dim>      dof_handler;
	
	SparsityPattern      sparsity_pattern;
	SparseMatrix<double> system_matrix;
	
	Vector<double>       solution;
	Vector<double>       system_rhs;
	
	ConvergenceTable convergence_table;
	
	int job_id;
	struct lapprm<dim> &lapp;
	
	const Point<dim> *grid_hrect_p1, *grid_hrect_p2; // \omega (which is a hyper rectangle)
	int basis_poly_deg; //degree of the tensor product polynomials (i.e. our basis functions) to be constructed
};

template <int dim>
void LaplaceSolve<dim>::setup_system()
{
	dof_handler.distribute_dofs(fe);
	
	sparsity_pattern.reinit(dof_handler.n_dofs(), dof_handler.n_dofs(),
	                        dof_handler.max_couplings_between_dofs());
	DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern);
	sparsity_pattern.compress();
	
	system_matrix.reinit(sparsity_pattern);
	
	solution.reinit(dof_handler.n_dofs());
	system_rhs.reinit(dof_handler.n_dofs());
}

template <int dim>
void LaplaceSolve<dim>::assemble_system()
{
	QGauss<dim> quadrature_formula((basis_poly_deg == 1) ? 2:3);
	
	const RightHandSide<dim> right_hand_side(lapp.Prob);
	
	FEValues<dim> fe_values(fe, quadrature_formula,
	                        update_values | update_gradients | update_quadrature_points | update_JxW_values);
	
	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points    = quadrature_formula.size();
	
	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	Vector<double>     cell_rhs(dofs_per_cell);
	
	std::vector<unsigned int> local_dof_indices(dofs_per_cell);
	
	typename DoFHandler<dim>::active_cell_iterator cell=dof_handler.begin_active(), endc=dof_handler.end();
	for (; cell!=endc; ++cell) {
		fe_values.reinit(cell);
		cell_matrix = 0;
		cell_rhs = 0;
		
		for (unsigned int q_point=0; q_point<n_q_points; ++q_point) {
			for (unsigned int i=0; i<dofs_per_cell; ++i) {
				for (unsigned int j=0; j<dofs_per_cell; ++j) {
					cell_matrix(i,j) += (fe_values.shape_grad(i, q_point) *
					                     fe_values.shape_grad(j, q_point) *
					                     fe_values.JxW(q_point));
				}
				
				double temp = fe_values.shape_value(i, q_point);
				Point<dim> p_temp = fe_values.quadrature_point(q_point);
				if (dim==2) {temp *= right_hand_side.value(p_temp[0], p_temp[1]);}
				if (dim==3) {temp *= right_hand_side.value(p_temp[0], p_temp[1], p_temp[2]);}
				temp *= fe_values.JxW(q_point);
				cell_rhs(i) += temp;
			}
		}
		
		cell->get_dof_indices(local_dof_indices);
		for (unsigned int i=0; i<dofs_per_cell; ++i) {
			for (unsigned int j=0; j<dofs_per_cell; ++j) {
				system_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i,j));
			}
			
			system_rhs(local_dof_indices[i]) += cell_rhs(i);
		}
	}
	
	std::map<unsigned int,double> boundary_values;
	VectorTools::interpolate_boundary_values(dof_handler, 0, BoundaryValues<dim>(job_id, lapp), boundary_values);
	MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution, system_rhs);
}

template <int dim>
void LaplaceSolve<dim>::solve()
{
	SolverControl solver_control(
		lapp.sp->subp[(lapp.sp->selective_sol == true) ? job_id:0].max_nof_iters,
		lapp.sp->subp[(lapp.sp->selective_sol == true) ? job_id:0].tol);
	SolverCG<> cg(solver_control);
	cg.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
	
	std::cout<<"   "<<solver_control.last_step()<<" CG iterations needed to obtain convergence."<<std::endl;
}

#define format__data_out(arg) \
do { \
	outfile << output_dir << "solution-" << dim << "d_" << job_id << "."#arg; \
	std::ofstream output(outfile.str().data()); \
	data_out.write_ ## arg (output); \
} while (0);

template <int dim>
void LaplaceSolve<dim>::output_results() const
{
	DataOut<dim> data_out;
	
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, "solution");
	
	data_out.build_patches();
	
	char output_dir[] = "../solution/";
	mkdir(output_dir, 0777);
	
	/* use gnuplot. e.g.:
	 * $ gnuplot
	 * gnuplot> splot 'solution-2d_0.gpl'
	 */
	#define OUTPUT_TYPE 4  /* 1:dx -opendx-, 2:gnuplot, 3:povray, 4:tecplot */
	std::ostringstream outfile;
	switch (OUTPUT_TYPE) {
		case 1: format__data_out(dx); break;
		case 2: format__data_out(gnuplot); break;
		case 3: format__data_out(povray); break;
		case 4: format__data_out(tecplot); break;
	}
	
	if (OUTPUT_TYPE == 2) { //gnuplot
		std::ostringstream command;
		command << "echo '" << "#$dim#$";
		for (int i=0; i<dim; i++) {command << " " << lapp.Prob.D[i];}
		command << "' >> " << outfile.str();
		system(command.str().data());
	}
}

template <int dim>
void LaplaceSolve<dim>::process_solution(const unsigned int cycle)
{
	Vector<float> difference_per_cell (triangulation.n_active_cells());
	
	VectorTools::integrate_difference (dof_handler,
										solution,
										Solution<dim>(lapp.Prob),
										difference_per_cell,
										QGauss<dim>(3),
										VectorTools::L2_norm);
	const double L2_error = difference_per_cell.l2_norm();
	
	//    VectorTools::integrate_difference (dof_handler,
	//                                      solution,
	//                                      Solution<dim>(lapp.Prob),
	//                                      difference_per_cell,
	//                                      QGauss<dim>(3),
	//                                      VectorTools::H1_seminorm);
	//    const double H1_error = difference_per_cell.l2_norm();
	const double H1_error = 0;
	
	const QTrapez<1>     q_trapez;
	const QIterated<dim> q_iterated (q_trapez, 5);
	VectorTools::integrate_difference (dof_handler,
										solution,
										Solution<dim>(lapp.Prob),
										difference_per_cell,
										q_iterated,
										VectorTools::Linfty_norm);
	const double Linfty_error = difference_per_cell.linfty_norm();
	
	const unsigned int n_active_cells=triangulation.n_active_cells();
	const unsigned int n_dofs=dof_handler.n_dofs();
	
	convergence_table.add_value("cycle", cycle);
	convergence_table.add_value("cells", n_active_cells);
	convergence_table.add_value("dofs", n_dofs);
	convergence_table.add_value("L2", L2_error);
	convergence_table.add_value("H1", H1_error);
	convergence_table.add_value("Linfty", Linfty_error);
}

#define print_convergence_table() \
do { \
	std::cout << std::endl; \
	convergence_table.set_precision("L2", 3); \
	convergence_table.set_precision("H1", 3); \
	convergence_table.set_precision("Linfty", 3); \
  \
	convergence_table.set_scientific("L2", true); \
	convergence_table.set_scientific("H1", true); \
	convergence_table.set_scientific("Linfty", true); \
  \
	convergence_table.set_tex_caption("cells", "\\# cells"); \
	convergence_table.set_tex_caption("dofs", "\\# dofs"); \
	convergence_table.set_tex_caption("L2", "L^2-error"); \
	convergence_table.set_tex_caption("H1", "H^1-error"); \
	convergence_table.set_tex_caption("Linfty", "L^\\infty-error"); \
  \
	convergence_table.set_tex_format("cells", "r"); \
	convergence_table.set_tex_format("dofs", "r"); \
  \
	convergence_table.write_text(std::cout); \
  \
	std::cout << std::endl; \
} while (0)

template <int dim>
void LaplaceSolve<dim>::run()
{
	for (unsigned int cycle=0; cycle < lapp.sp->subp[(lapp.sp->selective_sol == true) ? job_id:0].refine_times; ++cycle) {
		if (cycle == 0) { //make grid
			GridGenerator::hyper_rectangle(triangulation, *grid_hrect_p1, *grid_hrect_p2);
			std::cout<<"GridGenerator: hyper_rectangle: "<<std::endl;
			std::cout<<"\tfrom: ";  for (int i=0; i<dim; i++) {std::cout<<(*grid_hrect_p1)[i]<<" ";}  std::cout<<std::endl;
			std::cout<<"\tto: ";    for (int i=0; i<dim; i++) {std::cout<<(*grid_hrect_p2)[i]<<" ";}  std::cout<<std::endl;
		}
		
		triangulation.refine_global(1);
		setup_system();
		
		assemble_system();
		solve();
		
		process_solution(cycle);
	}
	
	print_convergence_table();
	
	output_results(); //output results of last cycle only
}

template <int dim>
void LaplaceSolve_main(struct lapprm<dim> *lapp, int job_id)
{
	std::cout<<"*========================================*"<<std::endl
	         <<"* Solving problem in "<<dim<<" space dimensions. *"<<std::endl
	         <<"*========================================*"<<std::endl<<std::endl;
	
	int basis_poly_deg; //degree of the tensor product polynomials (i.e. our basis functions) to be constructed
	for (basis_poly_deg = 1;  basis_poly_deg <= 2;  basis_poly_deg++) {
		std::cout << "Solving with Q"<<basis_poly_deg<<" elements, global refinement" << std::endl
		          << "===========================================" << std::endl << std::endl;
		
		LaplaceSolve<dim> laplace_problem(basis_poly_deg, job_id, *lapp);
		laplace_problem.run();
	}
}
//===================================================

#endif
