#include <iostream>

#include "ViennaCLSolver.h"

extern "C"
void viennacl_solve(const int m, const int nnz, const double *vals, int *rowidx, int *colptr, double *b, int *stat)
{
/*	std::cout << "Starting GPU Sparse Solve" << std::endl; */

	ViennaCLSolver solver;
		
	// the matrix is square & symetric, so we can ignore CSC-CSR conversion
	solver.solve(m, nnz, vals, rowidx, colptr, b, b, 1);

/*	std::cout << "GPU Sparse Solve Done" << std::endl;     */
	*stat = 0;
}
