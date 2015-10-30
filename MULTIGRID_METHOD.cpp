#include "MULTIGRID_METHOD.h"
#include "PROJECTION_2D.h"

void MULTIGRID_METHOD::Initialize(FIELD_STRUCTURE_2D<T>& solution_input, FIELD_STRUCTURE_2D<T>& rhs_input, FIELD_STRUCTURE_2D<int>& bc_input, FIELD_STRUCTURE_2D<T>& interface_levelset_input, POISSON_SOLVER_TYPE coarsest_level_poisson_solver_type_input)
{
	base_grid = &(solution_input.grid);
	num_levels = base_grid->RecommendMaxMultigridLevel(minimum_x_res);

	DeleteMemory();

	finest_level = 0;
	coarsest_level = num_levels - 1;

	projection_type = FREE_SURFACE_WATER;

	// Initialize list pointers of multigrid fields
	solution_multigrid.Initialize(num_levels);
	solution_multigrid_tmp.Initialize(num_levels);
	rhs_multigrid.Initialize(num_levels);
	res_multigrid.Initialize(num_levels);
	bc_multigrid.Initialize(num_levels);
	true_solution_multigrid.Initialize(num_levels);
	interface_levelset_multigrid.Initialize(num_levels);

	// Initialize instances of multigrid fileds except 0th ones received as parameters
	solution_multigrid[0] = &solution_input;
	solution_multigrid_tmp[0] = new FIELD_STRUCTURE_2D<T>();
	rhs_multigrid[0] = &rhs_input;
	bc_multigrid[0] = &bc_input;
	res_multigrid[0] = new FIELD_STRUCTURE_2D<T>();
	true_solution_multigrid[0] = new FIELD_STRUCTURE_2D<T>();
	interface_levelset_multigrid[0] = &interface_levelset_input;

	for (int level = 1; level < num_levels; level++)
	{
		solution_multigrid[level] = new FIELD_STRUCTURE_2D<T>();
		solution_multigrid_tmp[level] = new FIELD_STRUCTURE_2D<T>();
		rhs_multigrid[level] = new FIELD_STRUCTURE_2D<T>();
		bc_multigrid[level] = new FIELD_STRUCTURE_2D<int>();
		res_multigrid[level] = new FIELD_STRUCTURE_2D<T>();
		true_solution_multigrid[level] = new FIELD_STRUCTURE_2D<T>();
		interface_levelset_multigrid[level] = new FIELD_STRUCTURE_2D<T>();
	}

	// Initialize multi-grids
	int i_res = solution_input.grid.i_res, j_res = solution_input.grid.j_res;
	int i_start = solution_input.grid.i_start, j_start = solution_input.grid.j_start;
	T x_min = solution_input.grid.x_min, x_max = solution_input.grid.x_max;
	T y_min = solution_input.grid.y_min, y_max = solution_input.grid.y_max;

	for (int level = 1; level < num_levels; level++)
	{
		i_res /= 2; j_res /=2;
		solution_multigrid[level]->Initialize(i_res + 1, j_res + 1, i_start, j_start, x_min, y_min, x_max, y_max, 2, true, false, multithreading);
	}

	solution_multigrid_tmp[0]->Initialize(solution_multigrid[0]->grid, 2, multithreading);
	res_multigrid[0]->Initialize(solution_multigrid[0]->grid, 2, multithreading);
	true_solution_multigrid[0]->Initialize(solution_multigrid[0]->grid, 2, multithreading);

	for (int level = 1; level < num_levels; level++)
	{
		solution_multigrid_tmp[level]->Initialize(solution_multigrid[level]->grid, 2, multithreading);
		bc_multigrid[level]->Initialize(solution_multigrid[level]->grid, 2, multithreading);
		rhs_multigrid[level]->Initialize(solution_multigrid[level]->grid, 2, multithreading);
		res_multigrid[level]->Initialize(solution_multigrid[level]->grid, 2, multithreading);
		true_solution_multigrid[level]->Initialize(solution_multigrid[level]->grid, 2, multithreading);
		interface_levelset_multigrid[level]->Initialize(solution_multigrid[level]->grid, 2, multithreading);
	}

	// First time pressure initialization
	for (int level = 0; level < num_levels; level++)
	{
		solution_multigrid[level]->array_for_this.AssignAllValues((T)0);
	}

	poisson_solver_multigrid.Initialize(num_levels);

	T tolerance_level = tolerance;

	for (int level = 0; level < coarsest_level; level++)
	{
		poisson_solver_multigrid[level].Initialize(tolerance_level, max_iteration, &(solution_multigrid[level]->grid_ghost), 1, &(solution_multigrid[level]->partial_grids), num_levels, multithreading);
		poisson_solver_multigrid[level].gs_coefficients.Initialize(solution_multigrid[level]->grid, 2, multithreading);
		poisson_solver_multigrid[level].InitializeLinearSolver(CG);
	}

	poisson_solver_multigrid[coarsest_level].Initialize(tolerance_level, max_iteration, 0, 1, 0, num_levels, multithreading);
	poisson_solver_multigrid[coarsest_level].InitializeLinearSolver(CG);
}

void MULTIGRID_METHOD::InitializeFromBlock(const SCRIPT_BLOCK& multigrid_block)
{
	vcycle_iteration = multigrid_block.GetBoolean("vcycle_iteration", false);

	max_iteration = multigrid_block.GetInteger("max_iteration", 30);
	tolerance = multigrid_block.GetFloat("tolerance", (T)1e-4);
	sqr_tolerance = tolerance*tolerance;

	num_pre_smoothing = multigrid_block.GetInteger("num_pre_smoothing", 20);
	num_post_smoothing = multigrid_block.GetInteger("num_post_smoothing", 20);
	minimum_x_res = multigrid_block.GetInteger("minimum_x_res", 64);

	num_of_vcycle = multigrid_block.GetInteger("num_of_vcycle", 1);
}

void MULTIGRID_METHOD::Solve(const int& thread_id)
{
	PrepareForOneStep(thread_id);
	
	if (vcycle_iteration)
	{
		VCycleIteration(thread_id);
	}
	else
	{
		VCycle(thread_id);
	}
}

void MULTIGRID_METHOD::PrepareForOneStep(const int& thread_id)
{
	for (int level = 0; level <= coarsest_level; ++level)
	{
		if (test_number == 0)
		{
			BEGIN_GRID_ITERATION_2D(true_solution_multigrid[level]->partial_grids_ghost[thread_id])
			{
				T x_coor = true_solution_multigrid[level]->grid.x_min + i*true_solution_multigrid[level]->grid.dx, y_coor = true_solution_multigrid[level]->grid.y_min + j*true_solution_multigrid[level]->grid.dy;

				true_solution_multigrid[level]->array_for_this(i, j) = (POW2(x_coor) - POW2(POW2(x_coor)))*(POW2(POW2(y_coor)) - POW2(y_coor));
			}
			END_GRID_ITERATION_2D;
		}
	}

	for (int level = 0; level <= coarsest_level; ++level)
	{
		SetupBoundaryConditions(level, thread_id);
	}

	for (int level = 0; level < coarsest_level; ++level)
	{
		poisson_solver_multigrid[level].PrecomputeGSCoefficients(*(bc_multigrid[level]), thread_id);
	}
}

void MULTIGRID_METHOD::Smoothing(const int level, FIELD_STRUCTURE_2D<T>& solution, FIELD_STRUCTURE_2D<T>& solution_tmp, const FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& rhs, const int& max_itr, const int& thread_id)
{
	poisson_solver_multigrid[level].GaussSeidelSmoothing(solution, bc, rhs, max_itr, thread_id);
}

void MULTIGRID_METHOD::Restriction(const int& level, const int& thread_id)
{
	assert(level >= finest_level);
	assert(level < coarsest_level);

	poisson_solver_multigrid[level].CalculateResidual(*(res_multigrid[level]), *(solution_multigrid[level]), *(bc_multigrid[level]), *(rhs_multigrid[level]), thread_id);
	solution_multigrid[level + 1]->AssignAllValueGhost(thread_id, (T)0);

	ofstream fout;
	fout.open("residual_zero_level");
	for (int j = 0; j <= res_multigrid[0]->grid.j_end; j++)
	{
		for (int i = 0; i <= res_multigrid[0]->grid.i_end; i++)
		{
			fout << res_multigrid[0]->array_for_this(i, j) << " ";
		}
		fout << endl;
	}
	fout.close();

	//Injection(multithreading, *(res_multigrid[level]), *(rhs_multigrid[level + 1]), thread_id);
	FullWeighting(multithreading, *(res_multigrid[level]), *(rhs_multigrid[level + 1]), thread_id);
	
	// For interface levelset
	Injection(multithreading, *(interface_levelset_multigrid[level]), *(interface_levelset_multigrid[level + 1]), thread_id);

	fout.open("rhs_one_level");
	for (int j = 0; j <= res_multigrid[1]->grid.j_end; j++)
	{
		for (int i = 0; i <= res_multigrid[1]->grid.i_end; i++)
		{
			fout << rhs_multigrid[1]->array_for_this(i, j) << " ";
		}
		fout << endl;
	}
	fout.close();

	//Sampling(multithreading, *(res_multigrid[level]), *(rhs_multigrid[level + 1]), thread_id);
}

void MULTIGRID_METHOD::Prolongation(const int& level, const int& thread_id)
{
	assert(level >= finest_level);
	assert(level < coarsest_level);

	FillNonFullCellsWithZero(res_multigrid[level], bc_multigrid[level], thread_id);

	AdditiveSampling(multithreading, thread_id, *(solution_multigrid[level + 1]), *(solution_multigrid[level]), *(interface_levelset_multigrid[level]));

	ARRAY_2D<int> &bc_arr(bc_multigrid[level]->array_for_this);

	BEGIN_GRID_ITERATION_2D_ARR3(bc_multigrid[level]->partial_grids[thread_id], bc_multigrid[level]->array_for_this, solution_multigrid[level]->array_for_this, poisson_solver_multigrid[level].gs_coefficients.array_for_this)
	{
		if (bc_arr(arr1_ix) < 0)
		{
			(*(solution_multigrid[level]))(arr2_ix) = (T)0;
		}
		if (poisson_solver_multigrid[level].gs_coefficients(arr3_ix).center == (T)0)
		{
			(*(solution_multigrid[level]))(arr2_ix) = (T)0;
		}
	}
	END_GRID_ITERATION_2D;
}

void MULTIGRID_METHOD::SetupBoundaryConditions(const int& level, const int& thread_id)
{
	ARRAY_2D<int>			&bc_array(bc_multigrid[level]->array_for_this);
	GRID_STRUCTURE_2D		&grid(bc_multigrid[level]->grid);
	FIELD_STRUCTURE_2D<int> &bc(*(bc_multigrid[level]));
	FIELD_STRUCTURE_2D<T>	&solution(*(solution_multigrid[level]));
	FIELD_STRUCTURE_2D<T>   &true_solution(*(true_solution_multigrid[level]));

	BEGIN_GRID_ITERATION_2D(bc.partial_grids_ghost[thread_id])
	{
		if (i <= grid.i_start || i >= grid.i_end || j <= grid.j_start || j >= grid.j_end)
		{
			bc_array(i, j) = BC_DIR;
			solution(i, j) = true_solution(i, j);
			/*T x_coor = grid.x_min + i*grid.dx, y_coor = grid.y_min + j*grid.dy;

			solution(i, j) = (POW2(x_coor) - POW2(POW2(x_coor)))*(POW2(POW2(y_coor)) - POW2(y_coor));*/
		}
		else
		{
			bc_array(i, j) = BC_FULL;
		}
	}
	END_GRID_ITERATION_2D;
}

void MULTIGRID_METHOD::VCycleUpward(const int& thread_id)
{
	for (int level = coarsest_level - 1; level >= 0; level--)
	{
		Prolongation(level, thread_id);
		Smoothing(level, *(solution_multigrid[level]), *(solution_multigrid_tmp[level]), *(bc_multigrid[level]), *(rhs_multigrid[level]), num_post_smoothing, thread_id);
	}
	multithreading->Sync(thread_id);
	ofstream fout;
	fout.open("zero_level_four");
	for (int j = solution_multigrid[0]->grid.j_start - solution_multigrid[0]->ghost_width; j <= solution_multigrid[0]->grid.j_end + solution_multigrid[0]->ghost_width; j++)
	{
		for (int i = solution_multigrid[0]->grid.i_start - solution_multigrid[0]->ghost_width; i <= solution_multigrid[0]->grid.i_end + solution_multigrid[0]->ghost_width; i++)
		{
			fout << solution_multigrid[0]->array_for_this(i, j) << " ";
		}
		fout << endl; 
	}
	fout.close();
}

void MULTIGRID_METHOD::VCycle(const int& thread_id)
{
	VCycleDownward(thread_id);
	VCycleCoarsest(thread_id);
	VCycleUpward(thread_id);
}

void MULTIGRID_METHOD::VCycleIteration(const int& thread_id)
{
	for (int i = 0; i < num_of_vcycle; i++)
	{
		VCycleDownward(thread_id);
		VCycleCoarsest(thread_id);
		VCycleUpward(thread_id);
	}
}
void MULTIGRID_METHOD::VCycleCoarsest(const int& thread_id)
{
	int coarsest_level_resolution(solution_multigrid[coarsest_level]->grid.i_res);
	
	if (coarsest_level_resolution == 3)
	{
		BEGIN_HEAD_THREAD_WORK
		{
			solution_multigrid[coarsest_level]->array_for_this(1, 1) = solution_multigrid[coarsest_level]->array_for_this(1, 0) + solution_multigrid[coarsest_level]->array_for_this(1, 2) + solution_multigrid[coarsest_level]->array_for_this(0, 1) + solution_multigrid[coarsest_level]->array_for_this(2, 1) - solution_multigrid[coarsest_level]->dx*solution_multigrid[coarsest_level]->dx*rhs_multigrid[coarsest_level]->array_for_this(1, 1);
		}
		END_HEAD_THREAD_WORK;
	}
	else
	{
		poisson_solver_multigrid[coarsest_level].Solve(*(solution_multigrid[coarsest_level]), *(bc_multigrid[coarsest_level]), *(rhs_multigrid[coarsest_level]), thread_id);
		multithreading->Sync(thread_id);
	}
	
	ofstream fout;
	fout.open("coarsest_one");
	for (int j = solution_multigrid[coarsest_level]->grid.j_start; j <= solution_multigrid[coarsest_level]->grid.j_end; j++)
	{
		for (int i = solution_multigrid[coarsest_level]->grid.i_start; i <= solution_multigrid[coarsest_level]->grid.i_end; i++)
		{
			fout << solution_multigrid[coarsest_level]->array_for_this(i, j) << " ";
		}
		fout << endl; 
	}
	fout.close();
}

void MULTIGRID_METHOD::VCycleDownward(const int& thread_id)
{
	for (int level = 0; level < coarsest_level; level++)
	{
		Smoothing(level, *(solution_multigrid[level]), *(solution_multigrid_tmp[level]), *(bc_multigrid[level]), *(rhs_multigrid[level]), num_pre_smoothing, thread_id);
		Restriction(level, thread_id);
	}
	multithreading->Sync(thread_id);
	
	ofstream fout;
	fout.open("zero_level_one");
	for (int j = solution_multigrid[0]->grid.j_start; j <= solution_multigrid[0]->grid.j_end; j++)
	{
		for (int i = solution_multigrid[0]->grid.i_start; i <= solution_multigrid[0]->grid.i_end; i++)
		{
			fout << solution_multigrid[0]->array_for_this(i, j) << " ";
		}
		fout << endl; 
	}
	fout.close();
}

void MULTIGRID_METHOD::FillNonFullCellsWithZero(FIELD_STRUCTURE_2D<T>* scalar_field_input, FIELD_STRUCTURE_2D<int>* bc_field_input, const int& thread_id)
{
	GRID_STRUCTURE_2D	&grid(bc_field_input->grid);
	ARRAY_2D<int>		&bc_arr(bc_field_input->array_for_this);
	ARRAY_2D<T>			&scalar_arr(scalar_field_input->array_for_this);

	BEGIN_GRID_ITERATION_2D(bc_field_input->partial_grids_ghost[thread_id])
	{
		assert(bc_arr.i_end == scalar_arr.i_end);

		if (bc_arr(i, j) >= BC_FULL)
		{
			continue;
		}
		else
		{
			scalar_arr(i, j) = (T)0;
		}
	}
	END_GRID_ITERATION_2D;
}

void MULTIGRID_METHOD::DeleteMemory()
{
	if (solution_multigrid.length == 0)
	{
		return;
	}

	DELETE_POINTER(res_multigrid[0]);
	DELETE_POINTER(solution_multigrid_tmp[0]);
	
	for (int i = 1; i < num_levels; i++)
	{
		DELETE_POINTER(solution_multigrid[i]);
		DELETE_POINTER(rhs_multigrid[i]);
		DELETE_POINTER(res_multigrid[i]);
		DELETE_POINTER(bc_multigrid[i]);
		DELETE_POINTER(solution_multigrid_tmp[i]);
		DELETE_POINTER(true_solution_multigrid[i]);
	}
}