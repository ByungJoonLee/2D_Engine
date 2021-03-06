#pragma once

#include "LINEAR_SOLVER.h"
#include "LEVELSET_1D.h"
#include "LEVELSET_2D.h"
#include "CG_METHOD.h"
#include "PCG_METHOD.h"
#include "BICGSTAB_METHOD.h"
#include "GS_COEFFICIENTS.h"

class POISSON_SOLVER
{
public: // Essential Data
	T									tolerance;
	T									sqr_tolerance;
	int									max_iteration;
	
	int									ghost_width;

	GRID_STRUCTURE_2D*					grid_ghost;
	ARRAY<GRID_STRUCTURE_2D>*			partial_grids;

public: // Temporary
	int									num_levels;

public: // For matrix-vector type solver
	CSR_MATRIX<T>						A;
	VECTOR_ND<T>						x;
	VECTOR_ND<T>						b;

public: // pre-computed matrix-free Gauss-Seidel multiplication coefficients for direct solvers
	FIELD_STRUCTURE_2D<GS_COEFFICIENTS>	gs_coefficients;

public: // Solver
	LINEAR_SOLVER*						linear_solver;
	int									solver_type;

public: // Density 
	T									density_p, density_m;
	FIELD_STRUCTURE_1D<T>				one_over_density_1d;
	FIELD_STRUCTURE_2D<T>				one_over_density;

public: // Multithreading
	MULTITHREADING*						multithreading;

public: // Constructors and Destructor
	POISSON_SOLVER(void)
		: tolerance((T)1e-4), sqr_tolerance(tolerance*tolerance), max_iteration(100), linear_solver(0), multithreading(0), solver_type(CG)
	{}

	~POISSON_SOLVER(void)
	{
		DELETE_POINTER(linear_solver);
	}

public: // Initialization Functions
	void Initialize(const T& tolerance_input, const int& max_itr_input, const int ghost_width_input, MULTITHREADING* multithreading_input);
	void Initialize(const T& tolerance_input, const int& max_itr_input, GRID_STRUCTURE_2D* grid_ghost_input, const int ghost_width_input, ARRAY<GRID_STRUCTURE_2D>* partial_grids_input, int num_levels_input, MULTITHREADING* multithreading_input);
	void InitializeLinearSolver(const POISSON_SOLVER_TYPE linear_solver_type);
	
	void InitializePressure(FIELD_STRUCTURE_1D<T>& pressure);
	void InitializePressure(FIELD_STRUCTURE_2D<T>& pressure);
	void InitializePressure(FIELD_STRUCTURE_1D<T>& pressure, const int& thread_id);
	void InitializePressure(FIELD_STRUCTURE_2D<T>& pressure, const int& thread_id);

public: // Solver
	void Solve(FIELD_STRUCTURE_1D<T>& pressure, FIELD_STRUCTURE_1D<T>& density, FIELD_STRUCTURE_1D<int>& bc, const FIELD_STRUCTURE_1D<T>& div, const FIELD_STRUCTURE_1D<T>& variable, const LEVELSET_1D& levelset, const FIELD_STRUCTURE_1D<T>& jc_on_solution, const FIELD_STRUCTURE_1D<T>& jc_on_derivative, const int& thread_id);
	void Solve(FIELD_STRUCTURE_2D<T>& solution, FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& b);
	void Solve(FIELD_STRUCTURE_2D<T>& solution, FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& b, const int& thread_id);
	void Solve(FIELD_STRUCTURE_2D<T>& pressure, FIELD_STRUCTURE_2D<T>& density, FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& div, const LEVELSET_2D& levelset, const FIELD_STRUCTURE_2D<T>& jc_on_solution, const FIELD_STRUCTURE_2D<T>& jc_on_derivative);
	void Solve(FIELD_STRUCTURE_2D<T>& pressure, FIELD_STRUCTURE_2D<T>& density, FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& div, const FIELD_STRUCTURE_2D<T>& variable, const LEVELSET_2D& levelset, const FIELD_STRUCTURE_2D<T>& jc_on_solution, const FIELD_STRUCTURE_2D<T>& jc_on_derivative, const int& thread_id);
	void SolveForViscosity(FIELD_STRUCTURE_2D<T>& velocity, const FIELD_STRUCTURE_2D<T>& coef_1, const FIELD_STRUCTURE_2D<T>& coef_2, const FIELD_STRUCTURE_2D<T>& coef_3, const FIELD_STRUCTURE_2D<T>& coef_4, const FIELD_STRUCTURE_2D<T>& coef_5, FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& explicit_term);
	void SolveForAxisymmetric(FIELD_STRUCTURE_2D<T>& pressure, const FIELD_STRUCTURE_2D<T>& coef_1, const FIELD_STRUCTURE_2D<T>& coef_2, const FIELD_STRUCTURE_2D<T>& coef_3, const FIELD_STRUCTURE_2D<T>& coef_4, FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& div);

public: // Member Functions
	void SetTolerance(const T& tolerance_input);
	
	void BuildLinearSystem(CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_2D<T>& solution, const FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& RHS);
	void BuildLinearSystem(CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_2D<T>& solution, const FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& RHS, const int& thread_id);
	void BuildLinearSystemNodeForAxiSymmetric(CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_2D<T>& coef_1, const FIELD_STRUCTURE_2D<T>& coef_2, const FIELD_STRUCTURE_2D<T>& coef_3, const FIELD_STRUCTURE_2D<T>& coef_4, const FIELD_STRUCTURE_2D<T>& pressure, const FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& div);
	void BuildLinearSystemNodeForSemiImplicitViscosity(CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_2D<T>& coef_1, const FIELD_STRUCTURE_2D<T>& coef_2, const FIELD_STRUCTURE_2D<T>& coef_3, const FIELD_STRUCTURE_2D<T>& coef_4, const FIELD_STRUCTURE_2D<T>& coef_5, const FIELD_STRUCTURE_2D<T>& velocity, const FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& explicit_term);
	void BuildLinearSystemNodeJumpConditionVaribleCoefficient(CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_2D<T>& pressure, const FIELD_STRUCTURE_2D<T>& variable, const FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& div, const LEVELSET_2D& levelset, const FIELD_STRUCTURE_2D<T>& jc_on_solution, const FIELD_STRUCTURE_2D<T>& jc_on_derivative);
	void BuildLinearSystemNodeJumpConditionVaribleCoefficient(CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_1D<T>& pressure, const FIELD_STRUCTURE_1D<T>& variable, const FIELD_STRUCTURE_1D<int>& bc, const FIELD_STRUCTURE_1D<T>& div, const LEVELSET_1D& levelset, const FIELD_STRUCTURE_1D<T>& jc_on_solution, const FIELD_STRUCTURE_1D<T>& jc_on_derivative, const int& thread_id);
	void BuildLinearSystemNodeJumpConditionVaribleCoefficient(CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_2D<T>& pressure, const FIELD_STRUCTURE_2D<T>& variable, const FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& div, const LEVELSET_2D& levelset, const FIELD_STRUCTURE_2D<T>& jc_on_solution, const FIELD_STRUCTURE_2D<T>& jc_on_derivative, const int& thread_id);
	
	// This is amazing technique:)
	int AssignSequentialindicesToFullCells(const FIELD_STRUCTURE_1D<int>& bc, const int& thread_id);
	int AssignSequentialindicesToFullCells(const FIELD_STRUCTURE_2D<int>& bc);
	int AssignSequentialindicesToFullCells(const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id);
	
	int CountNonZeroElements(const FIELD_STRUCTURE_1D<int>& bc, const int& thread_id);
	int CountNonZeroElements(const FIELD_STRUCTURE_2D<int>& bc);
	int CountNonZeroElements(const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id);

	void VectorToGrid(const VECTOR_ND<T>& x, FIELD_STRUCTURE_1D<T>& pressure, const FIELD_STRUCTURE_1D<int>& bc, const int& thread_id);
	void VectorToGrid(const VECTOR_ND<T>& x, FIELD_STRUCTURE_2D<T>& pressure, const FIELD_STRUCTURE_2D<int>& bc);
	void VectorToGrid(const VECTOR_ND<T>& x, FIELD_STRUCTURE_2D<T>& pressure, const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id);

	// For multigrid
	void PrecomputeGSCoefficients(const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id); 
	void GaussSeidelSmoothing(FIELD_STRUCTURE_2D<T>& solution, const FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& rhs, const int& max_itr, const int& thread_id);
	void CalculateResidual(FIELD_STRUCTURE_2D<T>& residual, const FIELD_STRUCTURE_2D<T>& solution, const FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& rhs, const int& thread_id);
};