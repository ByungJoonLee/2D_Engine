#pragma once

#include "MULTILEVEL_SOLVER.h"

class MULTIGRID_METHOD : public MULTILEVEL_SOLVER
{
public: 
	int num_levels;
	T tolerance, sqr_tolerance;
	int max_iteration, minimum_x_res;

public: // Main Members
	ARRAY<FIELD_STRUCTURE_2D<T>*>		solution_multigrid;
	ARRAY<FIELD_STRUCTURE_2D<T>*>		solution_multigrid_tmp;
	ARRAY<FIELD_STRUCTURE_2D<T>*>		rhs_multigrid;
	
	ARRAY<FIELD_STRUCTURE_2D<T>*>		res_multigrid;
	ARRAY<FIELD_STRUCTURE_2D<int>*>		bc_multigrid;

	// For Interface Condition
	ARRAY<FIELD_STRUCTURE_2D<T>*>		interface_levelset_multigrid;

	// For reference solution
	ARRAY<FIELD_STRUCTURE_2D<T>*>		true_solution_multigrid;

	GRID_STRUCTURE_2D*					base_grid;

	ARRAY<POISSON_SOLVER>				poisson_solver_multigrid;

	MULTITHREADING*						multithreading;

public: // Options 
	bool								vcycle_iteration;
	
	// For Numerical Test
	int									test_number;

public: // Constants
	int									finest_level, coarsest_level;
	int									num_pre_smoothing, num_post_smoothing;
	int									num_of_vcycle;

public: // Constructor and Destructor
	MULTIGRID_METHOD(MULTITHREADING* multithreading_input)
		: multithreading(multithreading_input), solution_multigrid(0), solution_multigrid_tmp(0), rhs_multigrid(0), res_multigrid(0), bc_multigrid(0), true_solution_multigrid(0), interface_levelset_multigrid(0), 
		tolerance((T)1e-4), num_pre_smoothing(10), vcycle_iteration(false), num_of_vcycle(1), MULTILEVEL_SOLVER(multithreading_input), test_number(0)
	{}

	~MULTIGRID_METHOD()
	{
		DeleteMemory();
	}

public: // Initialization Function
	void Initialize(FIELD_STRUCTURE_2D<T>& solution_input, FIELD_STRUCTURE_2D<T>& rhs_input, FIELD_STRUCTURE_2D<int>& bc_input, FIELD_STRUCTURE_2D<T>& interface_level_input, POISSON_SOLVER_TYPE coarsest_level_poisson_solver_type_input);
	void InitializeFromBlock(const SCRIPT_BLOCK& outer_block);

public: // Solver
	void Solve(const int& thread_id);

public: // Functions
	void PrepareForOneStep(const int& thread_id);
	void Smoothing(const int level, FIELD_STRUCTURE_2D<T>& solution, FIELD_STRUCTURE_2D<T>& solution_tmp, const FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& rhs, const int& max_itr, const int& thread_id);
	void Restriction(const int& level, const int& thread_id);
	void Prolongation(const int& level, const int& thread_id);
	
	void FillNonFullCellsWithZero(FIELD_STRUCTURE_2D<T>* scalar_field_input, FIELD_STRUCTURE_2D<int>* bc_field_input, const int& thread_id);
	void VCycle(const int& thread_id);
	void VCycleIteration(const int& thread_id);
	void VCycleUpward(const int& thread_id);
	void VCycleCoarsest(const int& thread_id);
	void VCycleDownward(const int& thread_id);
	void SetupBoundaryConditions(const int& level, const int& thread_id);
	void DeleteMemory();

};