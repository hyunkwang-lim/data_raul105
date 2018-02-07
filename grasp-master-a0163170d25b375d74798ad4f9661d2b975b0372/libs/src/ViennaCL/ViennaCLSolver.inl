#pragma once

#include <iostream>
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/linalg/ilu.hpp"
#include "viennacl/linalg/ichol.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/bicgstab.hpp"

template<typename T>
inline void shift(T *x, const int shift, const int sz) {
	for (int i = 0; i < sz; ++i) 
		x[i] += shift;
}

inline void ViennaCLSolver::solve(
	const int m, const int nnz,
	const double	*vals,
	int		*cols,
	int		*rows,
	const double	*vec,
	double	*solution,
	const int firstIdx
	) 
{
	viennacl::compressed_matrix<double> a;
	viennacl::vector<double> b((size_t)m); 
	viennacl::vector<double> x((size_t)m);

	for (int i = 0; i < m; ++i) {
		b[i] = vec[i];
		x[i] = 0;
	}

	// viennacl exepctes indexing at 0
	shift(cols, -firstIdx, nnz); 
	shift(rows, -firstIdx, m+1);

	a.set(rows, cols, vals, (size_t)m, (size_t)m, (size_t)nnz);

	// recover indexing for rest of the code
	shift(cols, +firstIdx, nnz); 
	shift(rows, +firstIdx, m+1);

	// configure the matrix preconditioner to use
	viennacl::linalg::ilu0_precond< viennacl::compressed_matrix<double> >    pre(a, viennacl::linalg::ilu0_tag());
	// viennacl::linalg::ilut_precond< viennacl::compressed_matrix<double> >    pre(a, viennacl::linalg::ilut_tag(18, 1e-2));
	// viennacl::linalg::ichol0_precond< viennacl::compressed_matrix<double> >    pre(a, viennacl::linalg::ichol0_tag());

	// configure the iterative solver
	// viennacl::linalg::cg_tag solver(1e-16, 10000);
	viennacl::linalg::bicgstab_tag solver(1e-16, 6000, 3000);
	
	// do the actual solve
	// x = viennacl::linalg::solve(a, b, solver);
	x = viennacl::linalg::solve(a, b, solver, pre);
	
	//std::cout << "solver took " << solver.iters() << " iterations, error is " << solver.error() << std::endl;

	for (int i = 0; i < m; ++i) {
		solution[i] = x[i];
	}
}
