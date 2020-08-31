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
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/QR>

namespace nla3d {

namespace math {

    using namespace Eigen;
	typedef Eigen::SparseMatrix<double, Eigen::RowMajor> Eig_sparse;

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

//-------------------------------------------------------------------------------------------------------
 // relaxation steps before and after transmitting a problem to the next grid level
 class Relaxation {
 protected:
 	int n_relax_iters;
 public:
 	virtual void relaxation_steps(const Eig_sparse &A, const VectorXd &b, VectorXd &u) = 0;
 	virtual ~Relaxation() {};
 	Relaxation(int n_relax_iters);
 };

 class GaussSeidelRelaxation: public Relaxation {
 public:
 	virtual void relaxation_steps(const Eig_sparse &A, const VectorXd &b, VectorXd &u);
 	virtual ~GaussSeidelRelaxation() {};
 	GaussSeidelRelaxation(int n_relax_iters);
 };

 //------------------------------------------------------------------------------------------------------
 // the class represents the interface used by multigrid solver to compute the interpolation operator I
 // used to transmit the residual vector to the coarse grid as r_coarse = (I)^T * r
 // and to obtain error vector for current grid from the coarse grid one: e = I * e_coarse
 // and to compute the next level grid operator A_coarse = (I)^T * A * I
 class InterpolationOperator {
 protected:
 	//threshold for strength criterion
 	double theta;
 public:
 	virtual std::shared_ptr<Eig_sparse> get_operator(const Eig_sparse &A) = 0;
 	virtual ~InterpolationOperator() {};

 	InterpolationOperator(double theta);
 };

 class ClassicInterpolationOperator: public InterpolationOperator {
 private:
 	/*has one in row i and column j if node j has strong influence on node i
	 * otherwise zero*/
	Eigen::SparseMatrix<int> get_strong_influence_matrix(const Eig_sparse &A);
	/* divides all nodes in two subsets - fine and coarse. 
	 * nodes from the coarse subset remain in the next level grid while
	 * fine nodes are used on the current level only and will be computed 
	 * using coarse nodes values */
	void coarsen(Eigen::SparseMatrix<int> &S, std::vector<char> &Coarsening
		, std::vector<int> &coarse_grid_indexes);
 public:
 	virtual std::shared_ptr<Eig_sparse> get_operator(const Eig_sparse &A);
 	virtual ~ClassicInterpolationOperator() {};

 	ClassicInterpolationOperator(double theta);
 	ClassicInterpolationOperator();
 };

 class AggregationInterpolationOperator: public InterpolationOperator {
 private:
 	//for now it is assumed that with each call next level operator is returned. 
 	//the first call differs from the next ones as it corresponds to the finest level operator, 
 	// working with unknown vector containing three degrees of freedom for each node. 
 	// on the next levels each node contain six degrees of freedom
 	bool first_call;
 	// rigid body rotation and translation vectors computed from the dof coordinate vector
 	// passed into the constructor
 	Eigen::MatrixXd near_kernel;
 	/*has one in row i and column j if node j has strong influence on node i
	 * otherwise zero*/
	Eigen::SparseMatrix<int> get_strong_influence_matrix(const Eig_sparse &A, int node_dimension);
	/* divides all nodes in disjoint subsets called aggregates. each aggregate is 
	 * presented as node on the next level grid */
	void aggregate(std::vector<std::vector<int> > &aggregates
		, std::vector<bool> &undecided, const Eigen::SparseMatrix<int, RowMajor> &S);
 public:
 	/* each call modifies the near kernel vector, so in current implementation 
 	 * this method produces interpolation operator for i level only once.
 	 * each operator corresponds to corresponding grid level 
 	 * the operator itself described in paper: doi = 10.1145/2816795.2818081 */
 	virtual std::shared_ptr<Eig_sparse> get_operator(const Eig_sparse &A);
 	virtual ~AggregationInterpolationOperator() {};

 	AggregationInterpolationOperator(std::vector<double> &dof_coordinates, double theta=0.48);
 };

 //------------------------------------------------------------------------------------------------------

 class Preconditioner{
 public:
 	virtual ~Preconditioner() {};
 	virtual void precondition_step(const Eig_sparse &A, const VectorXd &b, VectorXd &x) = 0;
 };

 //------------------------------------------------------------------------------------------------------

 class MultigridSolver: public EquationSolver, public Preconditioner {
 private:
 	//precision
 	double eps;

 	//number of levels
 	//each level operations set is a recursive call 
 	int recursion_depth;
 	//number of iterations after which solving solving stops 
 	int iterations_limit;

 	//std::shared_ptr<Relaxation> relaxation;
 	Relaxation *relaxation;
 	//std::shared_ptr<InterpolationOperator> interpolation;
 	InterpolationOperator *interpolation;

 	/*interpolation operators needed
 	 * to transform the solution on coarser grid to finer scale*/
 	std::vector<std::shared_ptr<const Eig_sparse>> interpolation_operators;
 	/*matrix of coefficients A on different grid levels*/
 	std::vector<std::shared_ptr<const Eig_sparse>> coarse_equation_operators;

	void recursive_solver(const Eig_sparse &A, const VectorXd &b, VectorXd &u, int level);
	void direct_solve(MatrixXd &A, VectorXd &b, VectorXd &u);
	void solve(const Eig_sparse &A, const VectorXd &b, VectorXd &x);
 public:
 	int resulting_iterations;
 	double resulting_precision;

	 virtual ~MultigridSolver() {};
	 //these methods are inherited from EquationSolver interface. only solveEquations is used
	 virtual void solveEquations (math::SparseSymMatrix* matrix, double* rhs, double* unknowns);
	 virtual void factorizeEquations(math::SparseSymMatrix* matrix);
	 virtual void substituteEquations(math::SparseSymMatrix* matrix, double* rhs, double* unknowns);
	 virtual void precondition_step(const Eig_sparse &A, const VectorXd &b, VectorXd &x);
	 
	 void compute (const Eig_sparse &A);

	void set_relaxation_procedure(Relaxation *);
	void set_interpolation_procedure(InterpolationOperator *);

	MultigridSolver(int recursion_depth, int iterations_limit=-1, double eps = 1e-9) :
			 recursion_depth(recursion_depth), eps(eps), 
			 resulting_precision(0.), iterations_limit(iterations_limit), resulting_iterations(0) {};
 };
 
 class ConjugateGradientSolver : public EquationSolver {
 private:
	 double precision;
	 int iter_limit;
	 double resulting_precision;
	 int resulting_iterations;
	 
	 Preconditioner *preconditioner;

	 void solve(const Eig_sparse &A, const VectorXd &b, VectorXd &x);
 public:

	 ConjugateGradientSolver(double precision, int iter_limit):precision(precision), iter_limit(iter_limit), resulting_precision(0), resulting_iterations(0), preconditioner(nullptr) {};
	 
	 virtual ~ConjugateGradientSolver() {};
	 virtual void solveEquations (math::SparseSymMatrix* matrix, double* rhs, double* unknowns);
	 virtual void factorizeEquations(math::SparseSymMatrix* matrix);
 	 virtual void substituteEquations(math::SparseSymMatrix* matrix, double* rhs, double* unknowns);

	 double get_precision();

	 int get_iterations();

	 void set_preconditioner(Preconditioner *);
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
