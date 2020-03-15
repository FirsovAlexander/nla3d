// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "math/EquationSolver.h"

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
			iparm[2] = MKL_Get_Max_Threads();
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

		void MultigridSolver::solveEquations(math::SparseSymMatrix *matrix, double *rhs, double *unknowns) {
			std::vector<Triplet<double>> entries;
			Eig_sparse A(matrix->nRows(), matrix->nColumns());
			nEq = matrix->nColumns();

			for (uint16 i = 0; i < nEq; i++)
				for (uint16 j = 0; j < nEq; j++){
					double val = matrix->value(i + 1, j + 1);
					if(val != 0) entries.emplace_back(i, j, val);
				}

			A.setFromTriplets(entries.begin(), entries.end());
			VectorXd b(matrix->nColumns());
			VectorXd u = VectorXd::Zero(matrix->nColumns());

			for(int i = 0; i < matrix->nColumns(); i++)
				b[i] = rhs[i];

			u = this->solve(A, b, u);
			for(int i = 0; i < matrix->nColumns(); i++)
				unknowns[i] = u[i];
		}

		void MultigridSolver::factorizeEquations(math::SparseSymMatrix *matrix) {
			//no implementation yet
		}

		void MultigridSolver::substituteEquations(math::SparseSymMatrix *matrix, double *rhs, double *unknowns) {
			//no implementation yet
		}

		void MultigridSolver::recursive_solver(const Eig_sparse &A, const VectorXd &b, VectorXd &u,
												int level, int n_smoothing_iters) {
			std::cout<<level<<std::endl;
			static int level_index = 0;

			if (level == recursion_depth || A.cols() <= 500) {
				MatrixXd A_ = A.toDense();
				VectorXd b_ = b;
				direct_solve(A_, b_, u);
				return;
			}

			//std::cout << "rec" << std::endl;
			smooth(A, b, u, n_smoothing_iters);
			//std::cout << "smoothed" << std::endl;

			std::shared_ptr<Eig_sparse> interpol;
			std::shared_ptr<Eig_sparse> A_coarse;

			if (interpols.size() >= level) {
				A_coarse = coarse_equation_operators[level - 1];
				interpol = interpols[level - 1];
			} else {
				interpol = get_interpolation_operator(A);
				interpols.push_back(interpol);
				A_coarse = std::make_shared<Eig_sparse>((interpol->transpose() * A * (*interpol)));
				coarse_equation_operators.push_back(A_coarse);
			}

			VectorXd residual = interpol->transpose() * (b - A * u);
			VectorXd u_coarse = VectorXd::Zero(residual.size());

			//std::cout << "solving coarser" << std::endl;
			recursive_solver(*A_coarse, residual, u_coarse, level + 1, n_smoothing_iters);
			u += (*interpol) * u_coarse;
			smooth(A, b, u, n_smoothing_iters);

			if (levels.size() > 0 && level == levels[level_index]){
				residual = interpol->transpose() * (b - A * u);
				recursive_solver(*A_coarse, residual, u_coarse, level + 1, n_smoothing_iters);
				u += (*interpol) * u_coarse;
				smooth(A, b, u, n_smoothing_iters);

				level_index++;
				if(level_index == levels.size() - 1) level_index = 0;
			}

		}

//template<typename SmoothingSolver>
		std::shared_ptr<Eig_sparse>
		MultigridSolver::get_interpolation_operator(const Eig_sparse &A) {

			MatrixXi S = std::move(get_strong_influence_matrix(A));

			//  C/F coarsening
			std::vector<char> Coarsening(A.cols());// coarse and fine grid indexes defined after c/f coarsening
			std::vector<int> coarse_grid_indexes;
			int total_nodes = 0;

			while (total_nodes != Coarsening.size()) {
				int lambda = 0;
				int node_with_max_lambda = 0;

				for (int j = 0; j < A.outerSize(); j++) {
					if (Coarsening[j] != 0)
						continue;
					int S_and_U = 0, S_and_F = 0;
					for (Eig_sparse::InnerIterator it(A, j); it; ++it) {
						if (int row = it.row(); S(row, j)) {
							if (Coarsening[row] == 'F')
								S_and_F++;
							else if (Coarsening[row] == 0)
								S_and_U++;
						}
					}
					if (int l = S_and_U + 2 * S_and_F; l > lambda) {
						lambda = l;
						node_with_max_lambda = j;
					}
				}
				Coarsening[node_with_max_lambda] = 'C';
				coarse_grid_indexes.push_back(node_with_max_lambda);
				total_nodes++;
				for (Eig_sparse::InnerIterator it(A, node_with_max_lambda); it; ++it) {
					if (S(it.row(), node_with_max_lambda) && Coarsening[it.row()] == 0) {
						Coarsening[it.row()] = 'F';
						total_nodes++;
					}
				}
			}

			if (std::find(Coarsening.begin(), Coarsening.end(), 0) != Coarsening.end()) {
				std::cout << "no" << std::endl;
			}

			std::sort(coarse_grid_indexes.begin(), coarse_grid_indexes.end());
			std::map<int, int> fine_to_coarse;
			int coarse_order = 0;

			for (int &i : coarse_grid_indexes)
				fine_to_coarse[i] = coarse_order++;

			std::shared_ptr<Eig_sparse> I(new Eig_sparse(A.cols(), coarse_grid_indexes.size()));
			std::vector<double> positive_neighbours(A.cols()), positive_interpol_neighbours(A.cols()), negative_neighbours(
					A.cols()), negative_interpol_neighbours(A.cols());
			std::vector<Triplet<double >> entries;

			for (int j = 0; j < A.outerSize(); j++) {
				for (Eig_sparse::InnerIterator it(A, j); it; ++it) {
					if (int i = it.row(); Coarsening[i] == 'F') {
						if (i == j) continue;
						if (it.value() < 0) {
							negative_neighbours[i] += it.value();
							if (Coarsening[j] == 'C') {
								negative_interpol_neighbours[i] += it.value();
							}
						} else {
							positive_neighbours[i] += it.value();
							if (Coarsening[j] == 'C') {
								positive_interpol_neighbours[i] += it.value();
							}
						}
					}
				}
			}

			for (int j = 0; j < A.outerSize(); j++) {
				for (Eig_sparse::InnerIterator it(A, j); it; ++it) {
					if (Coarsening[j] == 'C') {
						if (int i = it.row(); i != j && Coarsening[i] == 'F')
							entries.emplace_back(i, fine_to_coarse[j], it.value() < 0 ? -negative_neighbours[i] /
																						negative_interpol_neighbours[i] *
																						it.value() / A.coeff(i, i) :
																	   -positive_neighbours[i] /
																	   positive_interpol_neighbours[i] * it.value() /
																	   A.coeff(i, i));
						else if (i == j) entries.emplace_back(i, fine_to_coarse[i], 1.);
					}
				}
			}

			I->setFromTriplets(entries.begin(), entries.end());
			I->makeCompressed();
			return I;
		}

//template<typename SmoothingSolver>
		VectorXd
		MultigridSolver::solve(const Eig_sparse &A, const VectorXd &b, const VectorXd &x,
								int n_smoothing_iters) {
			VectorXd u = x;
			double norm = (A * u - b).norm();
			std::cout << "residual norm before solving: " << norm << std::endl;
			//std::ofstream f("multigrid.csv");
			while (norm > eps) {
				this->resulting_iterations++;

				//recursive_solver(A, b, u, 1, n_smoothing_iters);
				smooth(A, b, u, n_smoothing_iters);
				norm = (A * u - b).norm();
				std::cout << "residual norm: " << norm << std::endl;
				//f << resulting_iterations << "\t" << norm << std::endl;

				if (resulting_iterations == iterations_limit) {
					std::cout << "Iterations limit exceeded" << std::endl;
					this->resulting_precision = norm;
					break;
				}
			}

			//f.close();

			this->resulting_precision = norm;

			return u;
		}

		MatrixXi MultigridSolver::get_strong_influence_matrix(const Eig_sparse &A) {
			MatrixXi S = MatrixXi::Zero(A.rows(), A.cols());
			for (int i = 0; i < A.rows(); i++) {
				double max_positive = 0;
				double min_negative = 0;
				for (int j = 0; j < A.cols(); j++) {
					if (i != j) {
						if (double d = A.coeff(i, j); d < min_negative)
							min_negative = d;
						else if (d > max_positive)
							max_positive = d;
					}
				}
				for (int j = 0; j < A.cols(); j++) {
					if (i != j) {
						if (double d = A.coeff(i, j); (min_negative != 0 && d / min_negative >= theta) ||
													  (max_positive != 0 && d / max_positive >= theta))
							S(i, j) = 1;
					}
				}
			}

			return S;
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

		void
		MultigridSolver::smooth(const Eig_sparse &A, const VectorXd &b, VectorXd &u, int n_smoothing_iters) {
			for (int i = n_smoothing_iters; i >= 0; i--) {
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
	} //namespace math

} //namespace nla3d
