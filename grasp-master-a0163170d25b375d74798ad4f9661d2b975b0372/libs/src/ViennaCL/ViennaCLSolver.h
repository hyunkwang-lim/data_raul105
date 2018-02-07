#pragma once

// \brief A sparse matrix solver utilizing http://viennacl.sourceforge.net/.
// \author Michael Aspetsberger
class ViennaCLSolver {
public:
	// \brief solves the linear system A*solution=vec
	// \param m the size of the square sparse matrix
	// \param nnz the amount of non-zero elements
	// \param vals array of size nnz, holding all non-zero values of A in CSR format
	// \param cols array of size nnz, holding the columns of A in CSR format
	// \param rows array of size m+1, holding the row-pointers of A in CSR format, where the last entry is nnz+1
	// \param vec array of size m, the right hand side vector
	// \param solution output array of size m, the computed solution
	// \param firstIdx defines where indexing starts, e.g. 1 for fortran
	inline void solve(
		const int m, const int nnz,
		const double	*vals,
		int		*cols,
		int		*rows,
		const double	*vec,
		double	*solution,
		const int firstIdx = 1
		);
};

// inlining the implemenation gives a significant (~10%) boost with g++
#include "ViennaCLSolver.inl"