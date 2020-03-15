// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"
#include "math/Mat.h"
#include "math/SparseMatrix.h"
#include <Eigen/Core>
#include <memory>
#include <Eigen/Sparse>

namespace nla3d {

namespace math {

    using namespace Eigen;
	typedef Eigen::SparseMatrix<double, RowMajor> Eig_sparse;

class SparseSymMatrix;

// EquationSolver - abstract class for solving a system of linear equations.
// This class primarly used in FESolver class.

class EquationSolver {
public:
  virtual ~EquationSolver() { };
  virtual void solveEquations(math::SparseSymMatrix* matrix, double* rhs, double* unknowns) = 0;
  virtual void factorizeEquations(math::SparseSymMatrix* matrix) = 0;
  virtual void substituteEquations(math::SparseSymMatrix* matrix, double* rhs, double* unknowns) = 0;
  void setSymmetric (bool symmetric = true);
  void setPositive (bool positive = true);
protected:
  uint32 nEq = 0;

  // number of rhs 
	int nrhs = 1; 
  bool isSymmetric = true;
  bool isPositive = true;
};

class GaussDenseEquationSolver : public EquationSolver {
public:
  virtual ~GaussDenseEquationSolver() { };
  virtual void solveEquations (math::SparseSymMatrix* matrix, double* rhs, double* unknowns);
  virtual void factorizeEquations(math::SparseSymMatrix* matrix);
  virtual void substituteEquations(math::SparseSymMatrix* matrix, double* rhs, double* unknowns);
  static bool _solve(double* X, double* A, double* B, int n);
protected:

  dMat matA = dMat(1, 1);
};

 class MultigridSolver: public EquationSolver {
 private:
 	double theta, eps;

 	int recursion_depth;
 	int iterations_limit;
 	int resulting_iterations = 0;
 	double resulting_precision = 0;


 	/*this is needed for w and f multigrid cycles like
 	 *          *
 	  \  *  *  /
 	   \/ \/ \/      */
 	std::vector<int> levels;//indexes of levels on which coarsening continues till recursion_depth

 	/*interpolation operators needed
 	 * to transform the solution on coarser grid to finer scale*/
 	std::vector<std::shared_ptr<Eig_sparse>> interpols;
 	std::vector<std::shared_ptr<Eig_sparse>> coarse_equation_operators;

	 VectorXd solve(const Eig_sparse &A, const VectorXd &b, const VectorXd &x, int n_smoothing_iters = 4);
	 void recursive_solver(const Eig_sparse &A, const VectorXd &b, VectorXd &u, int level, int n_smoothing_iters);
	 void smooth(const Eig_sparse &A, const VectorXd &b, VectorXd &u, int n_smoothing_iters);
	 Eigen::SparseMatrix<int> get_strong_influence_matrix(const Eig_sparse &A);
	 void coarsen(Eigen::SparseMatrix<int> &S, std::vector<char> &Coarsening, std::vector<int> &coarse_grid_indexes);
	 std::shared_ptr<Eig_sparse> get_interpolation_operator(const Eig_sparse &A);
	 void direct_solve(MatrixXd &A, VectorXd &b, VectorXd &u);
 public:
 	virtual ~MultigridSolver() {};
	 /*MultigridSolver(int recursion_depth, int iterations_limit, double theta = 0.25, double eps = 1e-9) :
			 recursion_depth(recursion_depth), theta(theta), eps(eps), iterations_limit(iterations_limit) {};*/
	 MultigridSolver(int recursion_depth, int iterations_limit, double theta = 0.25
	 		, double eps = 1e-9, std::vector<int> levels=std::vector<int>()) :
			 recursion_depth(recursion_depth), theta(theta), eps(eps)
			 , iterations_limit(iterations_limit), levels(std::move(levels)) {};
 	virtual void solveEquations (math::SparseSymMatrix* matrix, double* rhs, double* unknowns);
 	virtual void factorizeEquations(math::SparseSymMatrix* matrix);
 	virtual void substituteEquations(math::SparseSymMatrix* matrix, double* rhs, double* unknowns);
 };

#ifdef NLA3D_USE_MKL
class PARDISO_equationSolver : public EquationSolver {
public:
  virtual ~PARDISO_equationSolver();
  virtual void solveEquations (math::SparseSymMatrix* matrix, double* rhs, double* unknowns);
  virtual void factorizeEquations(math::SparseSymMatrix* matrix);
  virtual void substituteEquations(math::SparseSymMatrix* matrix, double* rhs, double* unknowns);
protected:
  void initializePARDISO (math::SparseSymMatrix* matrix);
  void releasePARDISO ();

  // Internal solver memory pointer pt
	void *pt[64]; 
  // Paramaters for PARDISO solver (see MKL manual for clarifications)
	int iparm[64];

  // maximum number of numerical factorizations
	int maxfct = 1; 

  // which factorization to use
	int mnum = 1; 
  // don't print statistical information in file
	int msglvl = 0; 
  
  bool firstRun = true;
  // real symmetric undifinite defined matrix
	int mtype = -2; 
};
#endif //NLA3D_USE_MKL

extern EquationSolver* defaultEquationSolver;

} // namespace math

} //namespace nla3d
