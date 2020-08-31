// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "math/EquationSolver.h"

#include <chrono>
#include <Eigen/LU>
#include <Eigen/SparseLU>
#include <Eigen/Eigenvalues>

#ifdef NLA3D_USE_MKL

#include <mkl.h>
#include <vector>

#endif //NLA3D_USE_MKL

namespace nla3d {

	namespace math {

		EquationSolver *defaultEquationSolver = new GaussDenseEquationSolver;

		void EquationSolver::setSymmetric(bool symmetric) {
			isSymmetric = symmetric;
		}

		void EquationSolver::setPositive(bool positive) {
			isPositive = positive;
		}


		void GaussDenseEquationSolver::solveEquations(math::SparseSymMatrix *matrix, double *rhs, double *unknowns) {
			TIMED_SCOPE(t, "solveEquations");
			factorizeEquations(matrix);
			substituteEquations(matrix, rhs, unknowns);
		}


		void GaussDenseEquationSolver::factorizeEquations(math::SparseSymMatrix *matrix) {
			// nothing to do here..
			nEq = matrix->nRows();
			CHECK(nEq < 1000) << "GaussDenseEquationSolver works only with number of equations less than 1000";
		}


		void GaussDenseEquationSolver::substituteEquations(math::SparseSymMatrix *matrix,
														   double *rhs, double *unknowns) {
			// NOTE: actually here we perform all steps factorization and substitution
			// because this is simplest solver dedicated to perform functional tests
			//
			CHECK(nrhs == 1) << "GaussDenseEquationSolver support only 1 set of rhs values";


			if (matA.dM() == nEq && matA.dN() == nEq) {
				matA.zero();
			} else {
				matA.resize(nEq, nEq);
			}

			// fill dense matA from sparse matrix
			for (uint16 i = 0; i < nEq; i++)
				for (uint16 j = 0; j < nEq; j++)
					// NOTE: SparseMatrix getters works with indexes started from 1
					// TODO: we can fill dense matrix from sparse one in a more optimized way
					matA[i][j] = matrix->value(i + 1, j + 1);
			bool res = _solve(unknowns, matA.ptr(), rhs, nEq);
			CHECK(res == true) << "ERROR during solution";
		}

		bool GaussDenseEquationSolver::_solve(double *X, double *A,
											  double *B, int n) {
			// Gaussian elimination, with partial pivoting. It's an error if the
			// matrix is singular, because that means two constraints are
			// equivalent.
			int i, j, ip, jp, imax;
			double max, temp;

			for (i = 0; i < n; i++) {
				// We are trying eliminate the term in column i, for rows i+1 and
				// greater. First, find a pivot (between rows i and N-1).
				max = 0;
				for (ip = i; ip < n; ip++) {
					//if(ffabs(A[ip][i]) > max) {
					if (fabs(A[ip * n + i]) > max) {
						imax = ip;
						//max = ffabs(A[ip][i]);
						max = fabs(A[ip * n + i]);
					}
				}
				// Don't give up on a singular matrix unless it's really bad; the
				// assumption code is responsible for identifying that condition,
				// so we're not responsible for reporting that error.
				if (fabs(max) < 1e-20) return false;

				// Swap row imax with row i
				for (jp = 0; jp < n; jp++) {
					//SWAP(double, A[i][jp], A[imax][jp]);
					double tmp;
					tmp = A[i * n + jp];
					A[i * n + jp] = A[imax * n + jp];
					A[imax * n + jp] = tmp;
				}
				//SWAP(double, B[i], B[imax]);
				double tmp;
				tmp = B[i];
				B[i] = B[imax];
				B[imax] = tmp;

				// For rows i+1 and greater, eliminate the term in column i.
				for (ip = i + 1; ip < n; ip++) {
					//temp = A[ip][i]/A[i][i];
					temp = A[ip * n + i] / A[i * n + i];

					for (jp = i; jp < n; jp++) {
						//A[ip][jp] -= temp*(A[i][jp]);
						A[ip * n + jp] -= temp * (A[i * n + jp]);
					}
					B[ip] -= temp * B[i];
				}
			}

			// We've put the matrix in upper triangular form, so at this point we
			// can solve by back-substitution.
			for (i = n - 1; i >= 0; i--) {
				//if(fabs(A[i][i]) < 1e-20) return false;
				if (fabs(A[i * n + i]) < 1e-20) return false;

				temp = B[i];
				for (j = n - 1; j > i; j--) {
					//temp -= X[j]*A[i][j];
					temp -= X[j] * A[i * n + j];
				}
				//X[i] = temp / A[i][i];
				X[i] = temp / A[i * n + i];
			}

			return true;
		}


#ifdef NLA3D_USE_MKL

		PARDISO_equationSolver::~PARDISO_equationSolver() {
			releasePARDISO();
		}

		void PARDISO_equationSolver::solveEquations(math::SparseSymMatrix *matrix, double *rhs, double *unknowns) {
			TIMED_SCOPE(t, "solveEquations");
			mkl_set_num_threads(0);

			factorizeEquations(matrix);
			substituteEquations(matrix, rhs, unknowns);
		}


		void PARDISO_equationSolver::substituteEquations(math::SparseSymMatrix *matrix,
														 double *rhs, double *unknowns) {

			//Back substitution and iterative refinement

			// initialize error code
			int error = 0;
			int phase = 33;
			int n = static_cast<int> (nEq);
			PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, matrix->getValuesArray(),
					(int *) matrix->getIofeirArray(),
					(int *) matrix->getColumnsArray(),
					NULL, &nrhs, iparm, &msglvl, rhs, unknowns, &error);

			CHECK(error == 0) << "ERROR during solution. Error code = " << error;
		}


		void PARDISO_equationSolver::factorizeEquations(math::SparseSymMatrix *matrix) {
			TIMED_SCOPE(t, "factorizeEquations");
			if (firstRun) {
				initializePARDISO(matrix);
			}
			firstRun = false;

			CHECK(nEq == matrix->nRows());

			int phase;

			// initialize error code
			int error = 0;

			// phase 22 is the numerical factorization
			phase = 22;
			int n = static_cast<int> (nEq);

			PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, matrix->getValuesArray(),
					(int *) matrix->getIofeirArray(),
					(int *) matrix->getColumnsArray(),
					NULL, &nrhs, iparm, &msglvl, NULL, NULL, &error);
			CHECK(error == 0) << "ERROR during numerical factorization. Error code = " << error;

		}

		void PARDISO_equationSolver::initializePARDISO(math::SparseSymMatrix *matrix) {
			for (uint16 i = 0; i < 64; i++) {
				iparm[i] = 0;
			}

			for (uint16 i = 0; i < 64; i++) {
				pt[i] = 0;
			}

			iparm[0] = 1; //no solver default
			iparm[1] = 2; //fill-in reordering from meris
			iparm[2] = 1;//MKL_Get_Max_Threads();
			iparm[3] = 0; //no iterative-direct algorithm
			iparm[4] = 0; //no user fill-in reducing permutation
			iparm[5] = 0; //write solution into x
			iparm[6] = 16; //default logical fortran unit number for output
			iparm[7] = 2; //max numbers of iterative refinement steps
			iparm[9] = 13; //pertrub the pivor elements with 1e-13
			iparm[10] = 1; //use nonsymmetric permutation  and scaling MPS
			iparm[13] = 0; //output: number of perturbed pivots
			iparm[17] = -1; //output: number of nonzeros in the factor LU
			iparm[18] = -1; //output: MFLOPS for LU factorization
			iparm[19] = 0; //output: number of CG Iterations

			LOG_IF(!isSymmetric, FATAL) << "For now PARDISO_equationSolver doesn't support non-symmetric matrices";
			if (isPositive) {
				mtype = 2;
				LOG(INFO) << "EquationSolver will use positive symmetric solver";
			} else {
				LOG(INFO) << "EquationSolver will use non-positive symmetric solver";
				mtype = -2;
			}

			nEq = matrix->nRows();
			int n = static_cast<int> (nEq);

			// initialize error code
			int error = 0;

			int phase = 11;

			PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, matrix->getValuesArray(),
					(int *) matrix->getIofeirArray(),
					(int *) matrix->getColumnsArray(),
					NULL, &nrhs, iparm, &msglvl, NULL, NULL, &error);
			CHECK(error == 0) << "ERROR during symbolic factorization. Error code = " << error;
			LOG(INFO) << "Number of nonzeros in factors = " << iparm[17] << ", number of factorization MFLOPS = "
					  << iparm[18];
		}

		void PARDISO_equationSolver::releasePARDISO() {
			int phase = -1;
			int n = static_cast<int> (nEq);

			// initialize error code
			int error = 0;
			//Termination and release of memory
			PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, NULL, NULL, NULL, NULL, &nrhs, iparm, &msglvl, NULL, NULL,
					&error);
			LOG_IF (error != 0, WARNING) << "ERROR during PARDISO termination. Error code = " << error;
		}

#endif //NLA3D_USE_MKL
		//all work is perfomed with eigen matrices now
		void MultigridSolver::solveEquations(math::SparseSymMatrix *matrix, double *rhs, double *unknowns) {
			std::vector<Triplet<double>> entries;
			Eigen::SparseMatrix<double, RowMajor> A(matrix->nRows(), matrix->nColumns());

			for (int i = 0; i < matrix->nRows(); i++){
				for (int j = matrix->getIofeirArray()[i] - 1; j < matrix->getIofeirArray()[i + 1] - 1; j++){
					if (matrix->getValuesArray()[j] == 0) continue;
					int j_ = matrix->getColumnsArray()[j] - 1;
					entries.emplace_back(i, j_, matrix->getValuesArray()[j]);
					if(i != j_) entries.emplace_back(j_, i, matrix->getValuesArray()[j]);
				}
			}

			A.setFromTriplets(entries.begin(), entries.end());
			A.makeCompressed();

			VectorXd b(matrix->nColumns());
			for(int i = 0; i < matrix->nColumns(); i++)
				b[i] = rhs[i];

			VectorXd u = VectorXd::Zero(matrix->nColumns());

			this->compute(A);
			this->solve(A, b, u);
			for(int i = 0; i < matrix->nColumns(); i++)
				unknowns[i] = u[i];
		}

		void MultigridSolver::factorizeEquations(math::SparseSymMatrix *matrix) {
			//no implementation yet
		}

		void MultigridSolver::substituteEquations(math::SparseSymMatrix *matrix, double *rhs, double *unknowns) {
			//no implementation yet
		}

		void MultigridSolver::compute(const Eig_sparse &A){

			interpolation_operators.clear();
			coarse_equation_operators.clear();

			auto before = std::chrono::steady_clock::now();

			//first element is computed from argument
			std::shared_ptr<Eig_sparse> I = interpolation->get_operator(A);
			std::shared_ptr<Eig_sparse> A_(new Eig_sparse(I->cols(), I->cols()));
			*A_ = I->transpose() * A * (*I);

			interpolation_operators.push_back(I);
			coarse_equation_operators.push_back(A_);
			
			//rec_depth includes fine level, that is number of coarse levels is rec_depth - 1
			//next elements are computed from latest vector element
			for (int i = 0; i < recursion_depth - 1; i++){
				auto I_ = interpolation->get_operator(*coarse_equation_operators[i]);
				std::shared_ptr<Eig_sparse> coarser_A(new Eig_sparse(I->cols(), I->cols()));
				*coarser_A = I_->transpose() * (*coarse_equation_operators[i]) * (*I_);

				coarse_equation_operators.push_back(coarser_A);
				interpolation_operators.push_back(I_);

				if(coarse_equation_operators[i+1]->cols() <= 20) break;
			}

			auto after = std::chrono::steady_clock::now();
			std::chrono::milliseconds t = std::chrono::duration_cast<std::chrono::milliseconds>(after - before);

			std::cout << "Prepared " << interpolation_operators.size() << " levels in " << t.count()/1000. << " seconds" << std::endl;
		}

		void MultigridSolver::recursive_solver(const Eig_sparse &A, const VectorXd &b, VectorXd &u,
												int level) {

			if (level-1 == interpolation_operators.size()){
				MatrixXd A_ = A.toDense();
				VectorXd b_ = b;
				direct_solve(A_, b_, u);
				return;
			}

			relaxation->relaxation_steps(A, b, u);

			std::shared_ptr<const Eig_sparse> interpol = interpolation_operators[level - 1];
			std::shared_ptr<const Eig_sparse> A_coarse = coarse_equation_operators[level - 1];

			VectorXd residual = interpol->transpose() * (b - A * u);
			VectorXd u_coarse = VectorXd::Zero(residual.size());

			recursive_solver(*A_coarse, residual, u_coarse, level + 1);
			u.noalias() += (*interpol) * u_coarse;

			relaxation->relaxation_steps(A, b, u);
		}

		void
		MultigridSolver::solve(const Eig_sparse &A, const VectorXd &b, VectorXd &u) {
			//VectorXd u = x;
			// std::ofstream out("plain.txt");
			double norm = (A * u - b).norm();
			// out << norm << std::endl;
			// std::cout << "Problem dimension:" << A.cols() << std::endl;
			std::cout << "residual norm before solving: " << norm << std::endl;

			int iters = 0;
			auto before = std::chrono::steady_clock::now();
			while (norm > eps) {
				iters++;

				recursive_solver(A, b, u, 1);
				norm = (A * u - b).norm();
				// out << norm << std::endl;
				std::cout << "residual norm: " << norm << std::endl;

				if (iters == iterations_limit) {
					std::cout << "Iterations limit exceeded" << std::endl;
					this->resulting_iterations = iters;
					this->resulting_precision = norm;
					break;
				}
			}
			auto after = std::chrono::steady_clock::now();

			this->resulting_precision = norm;
			this->resulting_iterations = iters;

			auto t = std::chrono::duration_cast<std::chrono::milliseconds>(after - before);
			std::cout << "Solved in " << this->resulting_iterations << " iterations and " << t.count()/1000. <<"seconds" << std::endl;
		}

		void MultigridSolver::direct_solve(MatrixXd &A, VectorXd &b, VectorXd &u) {
			//direct
			for (int i = 0; i < A.cols(); i++) {
				if (A(i, i) == 0) {
					int replace_eq_id = 0;
					for (int j = i + 1; j < A.rows(); j++) {
						if (A(j, i) != 0) {
							replace_eq_id = j;
							break;
						}
					}
					double exchange = 0;
					for (int j = i; j < A.cols(); j++) {
						exchange = A(i, j);
						A(i, j) = A(replace_eq_id, j);
						A(replace_eq_id, j) = exchange;
					}
					exchange = b[i];
					b[i] = b[replace_eq_id];
					b[replace_eq_id] = exchange;
				}
				for (int j = i + 1; j < A.rows(); j++) {
					if (A(j, i) == 0) continue;
					double coeff = A(j, i) / A(i, i);
					for (int k = i; k < A.cols(); k++) {
						A(j, k) -= A(i, k) * coeff;
					}
					b[j] -= b[i] * coeff;
				}
			}
			//backward
			for (int i = A.cols() - 1; i >= 0; i--) {
				for (int j = i + 1; j < A.cols(); j++) {
					b[i] -= A(i, j) * u[j];
				}
				u[i] = b[i] / A(i, i);
			}
		}

		void MultigridSolver::set_relaxation_procedure(Relaxation *relaxation){
			this->relaxation = relaxation;
		}

		void MultigridSolver::set_interpolation_procedure(InterpolationOperator *interpolation){
			this->interpolation = interpolation;
		}

		void MultigridSolver::precondition_step(const Eig_sparse &A, const VectorXd &b, VectorXd &x) {
			if (interpolation_operators.empty()) compute(A);
			recursive_solver(A, b, x, 1);
		}

		// InterpolationOperator --------------------------------------------------------------------------------------------------------

		InterpolationOperator::InterpolationOperator(double theta): theta(theta) {};

		// ClassicInterpolationOperator -------------------------------------------------------------------------------------------------

		ClassicInterpolationOperator::ClassicInterpolationOperator(): InterpolationOperator(0.25) {}

		ClassicInterpolationOperator::ClassicInterpolationOperator(double theta): InterpolationOperator(theta) {}

		Eigen::SparseMatrix<int> 
		ClassicInterpolationOperator::get_strong_influence_matrix(const Eig_sparse &A) {
			Eigen::SparseMatrix<int> S(A.rows(), A.cols());
			std::vector<Triplet<int >> entries;

			for (int i = 0; i < A.outerSize(); i++){
				double max_positive = 0;
				double min_negative = 0;
				for(Eig_sparse::InnerIterator it(A, i); it; ++it){
					if(i != it.col()){
						if (it.value() < min_negative)
							min_negative = it.value();
						else if (it.value() > max_positive)
							max_positive = it.value();
					}
				}
				for(Eig_sparse::InnerIterator it(A, i); it; ++it){
					if(i != it.col()){
						if ((min_negative != 0 && it.value() / min_negative >= theta) ||
						(max_positive != 0 && it.value() / max_positive >= theta))
							entries.emplace_back(i, it.col(), 1);
					}
				}
			}

			S.setFromTriplets(entries.begin(), entries.end());
			S.makeCompressed();
			return S;
		}

		std::shared_ptr<Eig_sparse>
		ClassicInterpolationOperator::get_operator(const Eig_sparse &A) {

			auto S = get_strong_influence_matrix(A);

			//  C/F coarsening
			// coarse and fine grid indexes defined after c/f coarsening
			std::vector<char> Coarsening(A.cols());

			std::vector<int> coarse_grid_indexes;

			coarsen(S, Coarsening, coarse_grid_indexes);

			std::sort(coarse_grid_indexes.begin(), coarse_grid_indexes.end());
			std::map<int, int> fine_to_coarse;
			int coarse_order = 0;
			for (int &i : coarse_grid_indexes)
				fine_to_coarse[i] = coarse_order++;

			std::shared_ptr<Eig_sparse> I(new Eig_sparse(A.cols(), coarse_grid_indexes.size()));
			std::vector<Triplet<double >> entries;
			for (int i : coarse_grid_indexes){
				entries.emplace_back(i, fine_to_coarse[i], 1.);
			}

			for (int i = 0; i < A.outerSize(); i++){
				if (Coarsening[i] == 'C') continue;

				double positive_neighbours = 0, positive_interpol_neighbours = 0
						, negative_neighbours = 0, negative_interpol_neighbours = 0;

				for (Eig_sparse::InnerIterator it(A, i); it; ++it){
					if (it.col() == i) continue;

					if (it.value() > 0){
						positive_neighbours += it.value();
						if (Coarsening[it.col()] == 'C')
							positive_interpol_neighbours += it.value();
					} else{
						negative_neighbours += it.value();
						if (Coarsening[it.col()] == 'C')
							negative_interpol_neighbours += it.value();
					}
				}
				for (Eig_sparse::InnerIterator it(A, i); it; ++it){
					if(it.col() != i && Coarsening[it.col()] == 'C'){
						double interpol_coeff = it.value() > 0 ? -positive_neighbours / positive_interpol_neighbours
								: -negative_neighbours / negative_interpol_neighbours;
						interpol_coeff *= it.value() / A.coeff(i, i);
						entries.emplace_back(i, fine_to_coarse[it.col()], interpol_coeff);
					}
				}
			}

			I->setFromTriplets(entries.begin(), entries.end());
			I->makeCompressed();
			return I;
		}

		void ClassicInterpolationOperator::coarsen(Eigen::SparseMatrix<int> &S
				, std::vector<char> &Coarsening, std::vector<int> &coarse_grid_indexes) {
			int total_nodes = 0;

			while (total_nodes != Coarsening.size()) {
				/* finding the node that has strong influence on the maximum
				 * size subsets of fine grid nodes and undecided nodes. 
				 * lambda weight is computed for each node on this purpose*/
				int lambda = -1;
				int node_with_max_lambda = -1;

				for (int j = 0; j < S.outerSize(); j++){
					if (Coarsening[j] != 0) continue;

					//subsets of nodes that are fine grid or undecided yet 
					//on which the current inspected node has strong influence
					int S_and_U = 0, S_and_F = 0;
					for(Eigen::SparseMatrix<int>::InnerIterator it(S, j); it; ++it){
						if (Coarsening[it.row()] == 'F')
							S_and_F++;
						else if (Coarsening[it.row()] == 0)
							S_and_U++;
					}
					if (int l = S_and_U + 2 * S_and_F; l > lambda) {
						//the node affecting the maximum amount of undecided and fine nodes is memorised
						lambda = l;
						node_with_max_lambda = j;
					}
				}
				//this node becomes coarse grid node
				Coarsening[node_with_max_lambda] = 'C';
				coarse_grid_indexes.push_back(node_with_max_lambda);
				total_nodes++;
				//and all the nodes on which this one has strong influence become fine nodes. 
				//they won't participate in next levels unknown vectors
				for (Eigen::SparseMatrix<int>::InnerIterator it(S, node_with_max_lambda); it; ++it) {
					if (Coarsening[it.row()] == 0) {
						Coarsening[it.row()] = 'F';
						total_nodes++;
					}
				}
			}
		}
		// AggregationInterpolationOperator ------------------------------------------------------------------------------------------------

		std::shared_ptr<Eig_sparse>
		AggregationInterpolationOperator::get_operator(const Eig_sparse &A){
			//node dimension. i. e. how many following vector elements correspond to one grid node
			int node_size;
			if (first_call) {
				node_size = 3;
				first_call = false;
			} else {
				node_size = 6;
			}		

			auto S = get_strong_influence_matrix(A, node_size);
			assert (S.cols() == S.rows());

			std::vector<std::vector<int> > aggregates;
			std::vector<bool> undecided(S.cols(), true);		

			aggregate(aggregates, undecided, S);
			if(std::find(undecided.begin(), undecided.end(), true) != undecided.end()){
				std::cout << "error";
				exit(0);
			}
			
			std::shared_ptr<Eig_sparse> I(new Eig_sparse(A.rows(), aggregates.size() * near_kernel.cols()));
			MatrixXd new_near_kernel(aggregates.size()*near_kernel.cols(), near_kernel.cols());
			std::vector<Triplet<double >> entries;		

			for(int j = 0; j < aggregates.size(); j ++){
				MatrixXd Q_i(aggregates[j].size() * node_size, near_kernel.cols());
				std::vector<Triplet<double >> Q_entries;
				std::sort(aggregates[j].begin(), aggregates[j].end());
				for(int i = 0; i < aggregates[j].size(); i++){
					for (int n_size = 0; n_size < node_size; n_size++){
						for(int k_size = 0; k_size < near_kernel.cols(); k_size++){
							Q_i(i * node_size + n_size, k_size) = near_kernel(aggregates[j][i] * node_size + n_size, k_size);
						}
					}
				}

				Eigen::HouseholderQR<MatrixXd> qr;
				qr.compute(Q_i);
				MatrixXd Q = static_cast<MatrixXd>(qr.householderQ()).block(0, 0, aggregates[j].size()*node_size, near_kernel.cols());		

				MatrixXd R = Q.transpose() * Q_i;
				new_near_kernel.block(j * near_kernel.cols(), 0, near_kernel.cols(), near_kernel.cols()) = R;

				for(int i = 0; i < aggregates[j].size(); i++){
					for (int n_size = 0; n_size < node_size; n_size++){
						for(int k_size = 0; k_size < near_kernel.cols(); k_size++){
							if(Q(i * node_size + n_size, k_size) == 0) continue;
							entries.emplace_back(aggregates[j][i] * node_size + n_size, j * near_kernel.cols() + k_size, Q(i * node_size + n_size, k_size));
						}
					}
				}
			}

			Eig_sparse I_(A.rows(), aggregates.size() * near_kernel.cols());
			I_.setFromTriplets(entries.begin(), entries.end());
			I_.makeCompressed();		

			entries.clear();

			auto inv = A.diagonal().asDiagonal().inverse();
			Eig_sparse a = inv * A;
			double w = 3./4. / a.blueNorm();
			for(int i = 0; i < A.cols(); i++){
				double a_ii = A.coeff(i, i);
				entries.emplace_back(i, i, 1 - w);
				for (Eig_sparse::InnerIterator it(A, i); it; ++it){
					if (it.row() == i) continue;
					entries.emplace_back(i, it.col(), -w * it.value()/a_ii);
				}
			}
			Eig_sparse Sm(A.rows(), A.cols());
			Sm.setFromTriplets(entries.begin(), entries.end());
			Sm.makeCompressed();

			*I = Sm * I_;
			near_kernel = new_near_kernel;
			return I;
		}

		Eigen::SparseMatrix<int> 
		AggregationInterpolationOperator::get_strong_influence_matrix(const Eig_sparse &A, int n_size){
			
			//strong influence is measured between blocks, not just matrix elements
			Eigen::SparseMatrix<int> S(A.rows()/n_size, A.cols()/n_size);
			std::vector<Triplet<int >> entries;

			//diagonal submatrices of power -1/2 are cached here. 
			std::vector<std::shared_ptr<MatrixXd> > blocks(A.cols()/n_size);

			//searching for strongly affecting blocks in row i (i. e. nodes having strrong influence on the i-th node)
			for (int i = 0; i < A.rows(); i += n_size){
				//norms of blocks a_ii^(-1/2) * a_ij * a_jj^(-1/2)
				std::vector<double> norms(A.cols()/n_size, 0.);
				double max_norm = 0;

				std::shared_ptr<MatrixXd> a_ii;
				if (blocks[i/n_size] == nullptr){
					//MatrixXd a = A.block(i * n_size, i * n_size, n_size, n_size).toDense().pow(-0.5);
					blocks[i/n_size] = std::make_shared<MatrixXd> (A.block(i, i, n_size, n_size).toDense().pow(-0.5));
				}
				a_ii = blocks[i/n_size];

				//iterating over nodes, participating in the i-th node values
				for (int j = 0; j < A.cols(); j += n_size){
					if(j == i) continue;

					auto ij_block = A.block(i, j, n_size, n_size).toDense();
					if (ij_block.isZero(0)) continue;

					std::shared_ptr<MatrixXd> a_jj;
					if ( blocks[j/n_size] == nullptr){
						blocks[j/n_size] = std::make_shared<MatrixXd>(A.block(j, j, n_size, n_size).toDense().pow(-0.5));
					}
					a_jj = blocks[j/n_size];

					norms[j/n_size] = ((*a_ii) * ij_block * (*a_jj)).operatorNorm();
					if (norms[j/n_size] > max_norm) max_norm = norms[j/n_size];
				}
				//all blocks with measure exceeding some threshold considered strongly affecting on node i
				for(int k = 0; k < norms.size(); k++){
					if(norms[k]/max_norm >= theta) 
						entries.emplace_back(i/n_size, k, 1);
				}
			}

			S.setFromTriplets(entries.begin(), entries.end());
			S.makeCompressed();
			return S;
		}

		void AggregationInterpolationOperator::aggregate(std::vector<std::vector<int> > &aggregates, std::vector<bool> &undecided, const Eigen::SparseMatrix<int, RowMajor> &S){
			for(int i = 0; i < undecided.size(); i++){
				if(undecided[i]){
					bool subset = true;
					for (Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(S, i);it ; ++it){
						if (!undecided[it.col()]) {
							subset = false;
							break;
						}
					}
					if (subset) {
						std::vector<int> c;
						c.push_back(i);
						undecided[i] = false;
						for (Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(S, i);it ; ++it){
							c.push_back(it.col());
							undecided[it.col()] = false;
						}
						aggregates.push_back(c);
					}
				}
			}
			
			for (auto & aggregate : aggregates){
				std::vector<int> to_add;
				for (int &i : aggregate){
					for (Eigen::SparseMatrix<int, RowMajor>::InnerIterator it(S, i);it ; ++it){
						if (undecided[it.col()]) {
							to_add.push_back(it.col());
							undecided[it.col()] = false;
						}
					}
				}
				aggregate.insert(aggregate.end(), to_add.begin(), to_add.end());
			}
			//std::cout << "end" << std::endl;
			std::vector<int> too_weak;
			for(int i = 0; i < undecided.size(); i++){
				if (undecided[i]) {
					too_weak.push_back(i);
					undecided[i] = false;
				}
			}
			if (!too_weak.empty()) aggregates.push_back(too_weak);
		}

		AggregationInterpolationOperator::AggregationInterpolationOperator(std::vector<double> &dof_coords, double theta): InterpolationOperator(theta), first_call(true) {
			//near_kernel_vectors - each col is a vector of rotation or translation
			//for grid nodes
			//format is: [x1, y1, z1, ... , x_n, y_n, z_n]
			near_kernel = MatrixXd::Zero(dof_coords.size(), 6);
			for (int i = 0; i < dof_coords.size(); i += 3){
				near_kernel(i, 0) = 1;
				near_kernel(i+1, 1) = 1;
				near_kernel(i+2, 2) = 1;

				//(0, -z, y)
				near_kernel(i+1, 3) = -dof_coords[i+2];
				near_kernel(i+2, 3) = dof_coords[i+1];
				//(z, 0, -x)
				near_kernel(i, 4) = dof_coords[i+2];
				near_kernel(i+2, 4) = -dof_coords[i];
				//(-y, x, 0)
				near_kernel(i, 5) = -dof_coords[i+1];
				near_kernel(i+1, 5) = dof_coords[i];
			}
		}
		// Relaxation -------------------------------------------------------------------------------------

		Relaxation::Relaxation(int n_relax_iters): n_relax_iters(n_relax_iters) {}

		// GaussSeidelRelaxation --------------------------------------------------------------------------

		GaussSeidelRelaxation::GaussSeidelRelaxation(int n_relax_iters): Relaxation(n_relax_iters) {}

		void GaussSeidelRelaxation::relaxation_steps(const Eig_sparse &A, const VectorXd &b, VectorXd &u){
			for (int i = n_relax_iters; i >= 0; i--) {
				for(int row = 0; row < A.outerSize(); row++){
					u[row] = b[row];
					double a_ii = 0;
					for(Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(A,row); it; ++it){
						if(it.col() == row) {
							a_ii = it.value();
							continue;
						}
						u[row] -= it.value() * u[it.col()];
					}
					u[row] /= a_ii;
				}
			}
		}

		//-------------------------------------------------------------------------------------------------

		void ConjugateGradientSolver::solveEquations(math::SparseSymMatrix *matrix, double *rhs, double *unknowns) {
			std::vector<Triplet<double>> entries;
			Eigen::SparseMatrix<double, RowMajor> A(matrix->nRows(), matrix->nColumns());

			for (int i = 0; i < matrix->nRows(); i++){
				for (int j = matrix->getIofeirArray()[i] - 1; j < matrix->getIofeirArray()[i + 1] - 1; j++){
					if (matrix->getValuesArray()[j] == 0) continue;
					int j_ = matrix->getColumnsArray()[j] - 1;
					entries.emplace_back(i, j_, matrix->getValuesArray()[j]);
					if(i != j_) entries.emplace_back(j_, i, matrix->getValuesArray()[j]);
				}
			}

			A.setFromTriplets(entries.begin(), entries.end());
			A.makeCompressed();

			VectorXd b(matrix->nColumns());
			VectorXd u = VectorXd::Zero(matrix->nColumns());

			for(int i = 0; i < matrix->nColumns(); i++)
				b[i] = rhs[i];

			this->solve(A, b, u);
			for(int i = 0; i < matrix->nColumns(); i++)
				unknowns[i] = u[i];
		}

		void ConjugateGradientSolver::solve(const Eig_sparse &A, const VectorXd &b, VectorXd &x) {
			
			VectorXd r = b - A*x;
			VectorXd p;

			if (preconditioner != nullptr){
				VectorXd z = VectorXd::Zero(x.size());
				preconditioner->precondition_step(A, r, z);
				p = z;
			} else {
				p = r;
			}

			double r_norm = r.dot(p);
			//std::ofstream f("conj.csv");
			auto before = std::chrono::steady_clock::now();
		
			while(1){

				this->resulting_iterations++;
		
				VectorXd p_by_a = A * p;
				double alpha = r_norm / p.dot(p_by_a);
				x.noalias() += alpha * p;
				r.noalias() -= alpha * p_by_a;
				
				double r_norm_new;
				if (preconditioner != nullptr){
					VectorXd z = VectorXd::Zero(x.size());
					preconditioner->precondition_step(A, r, z);

					r_norm_new = r.dot(z);
					double beta = r_norm_new / r_norm;
					p.noalias() = z + beta * p;
				} else {
					r_norm_new = r.squaredNorm();
					double beta = r_norm_new / r_norm;
					p.noalias() = r + beta * p;
				}
				
				r_norm = r_norm_new;
		
				double precision_ = r.norm();
				std::cout<<precision_<<std::endl;

				if (resulting_iterations == this->iter_limit){
					std::cout<<"Iterations limit exceeded"<<std::endl;
					this->resulting_precision = precision_;
					break;
				}
				if (precision_ < this->precision){
					this->resulting_precision = precision_;
					break;
				}
		
			}
			auto after = std::chrono::steady_clock::now();

			auto t = std::chrono::duration_cast<std::chrono::milliseconds>(after - before);
			std::cout << "Solved in " << this->resulting_iterations << " iterations and " << t.count()/1000. <<"seconds" << std::endl;
		}

		/*void ConjugateGradientSolver::solve(const Eig_sparse &A, const VectorXd &b, VectorXd &x) {
			std::cout << "dimenstion: " << A.cols() << std::endl;

			std::ofstream out("conj.txt");
			
			VectorXd r = b - A*x;
			out << r.norm() << std::endl;
			VectorXd p = r;

			double r_norm = r.dot(p);
			//std::ofstream f("conj.csv");
			auto before = std::chrono::steady_clock::now();
		
			while(1){

				this->resulting_iterations++;
		
				VectorXd p_by_a = A * p;
				double alpha = r_norm / p.dot(p_by_a);
				x.noalias() += alpha * p;
				r.noalias() -= alpha * p_by_a;
				out << r.norm() << std::endl;
				//out << r.norm() << std::endl;
				//z.setZero(); //= VectorXd::Zero(x.size());
				//precond.solve(A, r, z, 2);
				// amg.cycle(r, z);
		
				double r_norm_new = r.squaredNorm();
				double beta = r_norm_new / r_norm;
				p.noalias() = r + beta * p;

				r_norm = r_norm_new;
		
				double precision_ = r.norm();
				std::cout<<precision_<<std::endl;

				if (resulting_iterations == this->iter_limit){
					std::cout<<"Iterations limit exceeded"<<std::endl;
					this->resulting_precision = precision_;
					break;
				}
				if (precision_ < this->precision){
					this->resulting_precision = precision_;
					break;
				}
		
			}
			auto after = std::chrono::steady_clock::now();

			auto t = std::chrono::duration_cast<std::chrono::milliseconds>(after - before);
			std::cout << "Solved in " << this->resulting_iterations << " iterations and " << t.count()/1000. <<"seconds" << std::endl;
		}*/

		void ConjugateGradientSolver::factorizeEquations(math::SparseSymMatrix *matrix) {

		}

		void
		ConjugateGradientSolver::substituteEquations(math::SparseSymMatrix *matrix, double *rhs, double *unknowns) {

		}

		void ConjugateGradientSolver::set_preconditioner(Preconditioner *preconditioner) {
			this->preconditioner = preconditioner;
		}

		double ConjugateGradientSolver::get_precision(){
		 	return resulting_precision;
	 	}

	 	int ConjugateGradientSolver::get_iterations(){
		 	return resulting_iterations;
	 	}

	} //namespace math

} //namespace nla3d
