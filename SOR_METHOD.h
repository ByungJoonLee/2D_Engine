#pragma once

#include "LINEAR_SOLVER.h"
#include "CSR_MATRIX.h"
#include "VECTOR_ND.h"
#include "LINEAR_SOLVER.h"
#include "FIELD_STRUCTURE_2D.h"

class SOR_METHOD : public LINEAR_SOLVER
{
public: // Typedef
	typedef LINEAR_SOLVER BASE;

public: // Using Keyword
	using BASE::tolerance;
	using BASE::sqr_tolerance;
	using BASE::residual;
	using BASE::max_iteration;
	using BASE::num_iteration;
	using BASE::multithreading;
    using BASE::omega;

public: // Constructor and Destructor
	SOR_METHOD(void)
	{}

	~SOR_METHOD(void)
	{}

public: 
	void Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id);

	void SORMethod(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const int& thread_id);
	T    SORStep(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const int& thread_id);
};