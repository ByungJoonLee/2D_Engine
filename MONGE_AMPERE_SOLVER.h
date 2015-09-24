#pragma once

#include "COMMON_DEFINITION.h"
#include "WORLD_DISCRETIZATION_2D.h"
#include "LEVELSET_2D.h"
#include "CSR_MATRIX.h"
#include "LINEAR_SOLVER.h"
#include "GAUSS_SEIDEL_METHOD.h"
#include "SOR_METHOD.h"
#include "POISSON_SOLVER.h"
#include "BICGSTAB_METHOD.h"
#include "MATRIX_2X2.h"
#include "MATRIX_3X3.h"

class MONGE_AMPERE_SOLVER
{
public: // World Discretization
	WORLD_DISCRETIZATION_2D*			world_discretization;

public: // Main variables
	FIELD_STRUCTURE_2D<T>				u;
	FIELD_STRUCTURE_2D<VT>				grad_u;

	FIELD_STRUCTURE_2D<VT>				true_solution;
	
	FIELD_STRUCTURE_2D<T>				density_x;
	FIELD_STRUCTURE_2D<T>				density_y;
	
	ARRAY<FIELD_STRUCTURE_2D<T>>		MA_array;
	FIELD_STRUCTURE_2D<T>				MA;

	FIELD_STRUCTURE_2D<T>				partial_p;
	FIELD_STRUCTURE_2D<T>				partial_q;

	FIELD_STRUCTURE_2D<int>				bc;
	FIELD_STRUCTURE_2D<T>				RHS;
	FIELD_STRUCTURE_2D<T>				RHS_h;

	// For target normal vector and target boundary points
	int									n_y, N_point;

	// For boundary condition setup
	ARRAY<T>							H_n;
	int									max_lower, max_upper;

public: // Convenient variables and references
	GRID_STRUCTURE_2D					base_grid;

public: // For multithreading
	MULTITHREADING*						multithreading;

public: // For Newton Method
	CSR_MATRIX<T>						jacobian;
	VECTOR_ND<T>						x;
	VECTOR_ND<T>						b;

	VECTOR_ND<T>						Jb;
	int									num_iteration, max_iteration_for_sub_linear_solver, max_iteration_for_Newton_Method;
	T									delta;
	
	T									alpha;

	enum POISSON_SOLVER_TYPE            sub_linear_solver_type;
	LINEAR_SOLVER*						sub_linear_solver;
	
	enum POISSON_SOLVER_TYPE			poisson_solver_type;
	POISSON_SOLVER						poisson_solver;

	T									tolerance;

	bool								is_regularized;

	T									l2_norm, max_norm;

public: // Test
	int									test_number;

public: // Constructor and Destructor
	MONGE_AMPERE_SOLVER(void)
		: multithreading(0), sub_linear_solver(0), is_regularized(false), n_y(32), l2_norm(0), max_norm(0)
	{}

	~MONGE_AMPERE_SOLVER(void)
	{}

public: // Initialization Function
	void InitializeFromBlock(const SCRIPT_BLOCK& monge_ampere_eqn_block, MULTITHREADING* multithreading_input)
	{
		multithreading = multithreading_input;

		// Grid Initialization
		base_grid.Initialize(world_discretization->world_grid);

		// Initialize Fields
		u.Initialize(base_grid, 1, multithreading);
		grad_u.Initialize(base_grid, 1, multithreading);

		true_solution.Initialize(base_grid, 1, multithreading);
		
		density_x.Initialize(base_grid, 1, multithreading);
		density_y.Initialize(base_grid, 1, multithreading);
				
		bc.Initialize(base_grid, 1, multithreading);
		RHS.Initialize(base_grid, 1, multithreading);
		RHS_h.Initialize(base_grid, 1, multithreading);

		// We limited here as compact scheme
		MA_array.Initialize(2);
		MA.Initialize(base_grid, 1, multithreading);

		partial_p.Initialize(base_grid, 1, multithreading);
		partial_q.Initialize(base_grid, 1, multithreading);

		MA_array[0].Initialize(base_grid, 1, multithreading);
		MA_array[1].Initialize(base_grid, 1, multithreading);

		// For practical reason
		delta = base_grid.dx2;

		test_number = monge_ampere_eqn_block.GetInteger("test_number", (int)1);
		tolerance = monge_ampere_eqn_block.GetFloat("tolerance", (T)0.001);
		max_iteration_for_sub_linear_solver = monge_ampere_eqn_block.GetInteger("max_iteration_for_sub_linear_solver", (int)100);
		max_iteration_for_Newton_Method = monge_ampere_eqn_block.GetInteger("max_iteration_for_Newton_Method", (int)100);
		
		n_y = monge_ampere_eqn_block.GetInteger("n_y", (int)32);
		N_point = monge_ampere_eqn_block.GetInteger("N_point", (int)100);
		is_regularized = monge_ampere_eqn_block.GetBoolean("is_regularized", (bool)false);

		cout << "-------------------Monge-Ampere Solver Started-------------------" << endl;
		cout << "Test Number : #" << test_number << endl;
		cout << "n_y = " << n_y << endl;
		cout << "N_point = " << N_point << endl;

		// For sub solver
		const char* sub_linear_solver_type_input = monge_ampere_eqn_block.GetString("sub_linear_solver_type", "Null");
		
		if (!strcmp(sub_linear_solver_type_input, "GS"))
		{
			sub_linear_solver_type = GS;
		}
		else if (!strcmp(sub_linear_solver_type_input, "SOR"))
		{
			sub_linear_solver_type = SOR;
		}
		else if (!strcmp(sub_linear_solver_type_input, "BICG"))
		{
			sub_linear_solver_type = BICG;
		}
		else
		{
			sub_linear_solver_type = GS;
		}

		InitializeLinearSolver(sub_linear_solver_type);
		
		if (sub_linear_solver_type == SOR)
		{
			sub_linear_solver->omega = monge_ampere_eqn_block.GetFloat("omega", (T)1);
		}

		cout << "-------------------For Sublinear Solver-------------------" << endl;
		cout << "tolerance: " << tolerance << endl;
		cout << "max_iteration: " << max_iteration_for_sub_linear_solver << endl;
		
		if (sub_linear_solver_type == SOR)
		{
			cout << "Omega: " << sub_linear_solver->omega << endl;
		}

		switch (sub_linear_solver_type)
		{
		case NO_SOLVER:
			cout << "sublinear solver type: " << "NO SOLVER" << endl;
			break;

		case GS:
			cout << "sublinear solver type: " << "GS" << endl;
			break;
		
		case SOR:
			cout << "sublinear solver type: " << "SOR" << endl;
			break;

		case BICG:
			cout << "sublinear solver type: " << "BICG" << endl;
			break;

		default:
			break;
		}

		// For having initial value of Newton Method
		const char* poisson_solver_type_input = monge_ampere_eqn_block.GetString("poisson_solver_type", "Null");
		if (!strcmp(poisson_solver_type_input, "CG"))
		{
			poisson_solver_type = CG;
		}
		else
		{
			poisson_solver_type = CG;
		}
		if (!strcmp(poisson_solver_type_input, "PCG"))
		{
			poisson_solver_type = PCG;
		}
		else
		{
			poisson_solver_type = CG;
		}

		alpha = monge_ampere_eqn_block.GetFloat("alpha", (T)1);

		cout << "-------------------For Newton's Method-------------------" << endl;
		cout << "tolerance: " << tolerance << endl;
		cout << "max_iteration: " << max_iteration_for_Newton_Method << endl;
		
		if (is_regularized == true)
		{
			cout << "is_regularized: " << "true" << endl;
		}
		else
		{
			cout << "is_regularized: " << "false" << endl;
		}

		switch (poisson_solver_type)
		{
		case NO_SOLVER:
			cout << "poisson solver type: " << "NO SOLVER" << endl;
			break;

		case CG:
			cout << "poisson solver type: " << "CG" << endl;
			break;

		case PCG:
			cout << "poisson solver type: " << "PCG" << endl;
			break;

		default:
			break;
		}

		// Initialize Poisson solvers
		poisson_solver.Initialize(tolerance, max_iteration_for_sub_linear_solver, 0, multithreading);
		poisson_solver.InitializeLinearSolver(poisson_solver_type);
	}

	void InitializeLinearSolver(const POISSON_SOLVER_TYPE linear_solver_type);

public: // Member Function
	void SolveThreaded(const int& thread_id)
	{
		//SetupDensity(thread_id);
		SetupInitialForNewtonMethod(thread_id);
		SetupBoundaryCondition(bc, thread_id);
		NewtonMethod(jacobian, x, b, thread_id);
		ComputeGradient();
		ComputeErrorInL2Norm(thread_id);
		ComputeErrorInMaxNorm(thread_id);

		BEGIN_HEAD_THREAD_WORK
		{
			if (test_number == 1)
			{
				cout << "L2 norm: " << l2_norm << endl;
				cout << "Max norm: " << max_norm << endl;
			}
			if (test_number == 2)
			{
				cout << "L2 norm: " << l2_norm << endl;
				cout << "Max norm: " << max_norm << endl;
			}
			if (test_number == 3)
			{
				cout << "L2 norm: " << l2_norm << endl;
				cout << "Max norm: " << max_norm << endl;
			}
		}
		END_HEAD_THREAD_WORK;
	}
	
	void ComputeErrorInL2Norm(const int& thread_id)
	{
		BEGIN_HEAD_THREAD_WORK
		{
			l2_norm = 0;
		}
		END_HEAD_THREAD_WORK;
		
		BEGIN_GRID_ITERATION_2D(grad_u.partial_grids[thread_id])
		{
			if (test_number == 2)
			{
				T grd_x = grad_u(i, j).x, grd_y = grad_u(i, j).y;

				if (grad_u.draw_for_this(i, j) == true)
				{
					if (POW2(0.8*grd_x - 0.2*grd_y) + POW2(0.2*grd_x - 0.6*grd_y) < POW2(0.44))
					{
						grad_u.draw_for_this(i, j) = true;
					}
					else
					{
						grad_u.draw_for_this(i, j) = false;
					}
				}
			}

			if (test_number == 3)
			{
				T x_coor = grad_u.x_min + i*grad_u.dx, y_coor = grad_u.y_min + j*grad_u.dy;

				if (grad_u.draw_for_this(i, j) == true)
				{
					T domain_region_1 = POW2(x_coor + 0.1) + POW2(y_coor);
					T domain_region_2 = POW2(x_coor - 0.1) + POW2(y_coor);
				
					if (domain_region_1 < POW2(0.85) && x_coor < -0.1)
					{
						grad_u.draw_for_this(i, j) = true;
					}
					else if (domain_region_2 < POW2(0.85) && x_coor > 0.1)
					{
						grad_u.draw_for_this(i, j) = true;
					}
					else
					{
						grad_u.draw_for_this(i, j) = false;
					}
				}
			}
		}
		END_GRID_ITERATION_2D;

		BEGIN_GRID_ITERATION_2D(grad_u.partial_grids[thread_id])
		{
			if (grad_u.draw_for_this(i, j) == true)
			{
				l2_norm += sqrt(POW2(true_solution(i, j).x - grad_u(i, j).x) + POW2(true_solution(i, j).y - grad_u(i, j).y))*grad_u.grid.dx2;
			}
			else
			{
				continue;
			}
		}
		END_GRID_ITERATION_SUM(l2_norm);
	}
	
	void ComputeErrorInMaxNorm(const int& thread_id)
	{
		BEGIN_HEAD_THREAD_WORK
		{
			max_norm = 0;
		}
		END_HEAD_THREAD_WORK;
		
		BEGIN_GRID_ITERATION_2D(grad_u.partial_grids[thread_id])
		{
			if (test_number == 2)
			{
				T grd_x = grad_u(i, j).x, grd_y = grad_u(i, j).y;

				if (grad_u.draw_for_this(i, j) == true)
				{
					if (POW2(0.8*grd_x - 0.2*grd_y) + POW2(0.2*grd_x - 0.6*grd_y) < POW2(0.44))
					{
						grad_u.draw_for_this(i, j) = true;
					}
					else
					{
						grad_u.draw_for_this(i, j) = false;
					}
				}
			}

			if (test_number == 3)
			{
				T x_coor = grad_u.x_min + i*grad_u.dx, y_coor = grad_u.y_min + j*grad_u.dy;

				if (grad_u.draw_for_this(i, j) == true)
				{
					T domain_region_1 = POW2(x_coor + 0.1) + POW2(y_coor);
					T domain_region_2 = POW2(x_coor - 0.1) + POW2(y_coor);
				
					if (domain_region_1 < POW2(0.85) && x_coor < -0.1)
					{
						grad_u.draw_for_this(i, j) = true;
					}
					else if (domain_region_2 < POW2(0.85) && x_coor > 0.1)
					{
						grad_u.draw_for_this(i, j) = true;
					}
					else
					{
						grad_u.draw_for_this(i, j) = false;
					}
				}
			}
		}
		END_GRID_ITERATION_2D;

		BEGIN_GRID_ITERATION_2D(grad_u.partial_grids[thread_id])
		{
			if (grad_u.draw_for_this(i, j) == true)
			{
				max_norm = MAX(max_norm, sqrt(POW2(true_solution(i, j).x - grad_u(i, j).x) + POW2(true_solution(i, j).y - grad_u(i, j).y)));
			}
			else
			{
				continue;
			}
		}
		END_GRID_ITERATION_MAX_2D(max_norm);
	}

	void AssignStencil(const FIELD_STRUCTURE_2D<T>& rho_x, const FIELD_STRUCTURE_2D<T>& rho_y_1, const FIELD_STRUCTURE_2D<T>& rho_y_2, const int& thread_id)
	{
		// Speedup Variable
		T one_over_dx2 = u.grid.one_over_dx2, one_over_dy2 = u.grid.one_over_dy2;
		T one_over_2dx = u.grid.one_over_2dx, one_over_2dy = u.grid.one_over_2dy;

		int mid_i = (T)0.5*(u.grid.i_start + u.grid.i_end);
		int mid_j = (T)0.5*(u.grid.j_start + u.grid.j_end);

		BEGIN_GRID_ITERATION_2D(u.partial_grids[thread_id])
		{
			T dx, dy, dxx, dyy, dv, dvp, dvv, dvpvp;
			
			// x, y derivatives
			dx = one_over_2dx*(u(i + 1, j) - u(i - 1, j));
			dy = one_over_2dy*(u(i, j + 1) - u(i, j - 1));
			dxx = one_over_dx2*(u(i + 1, j) - 2*u(i, j) + u(i - 1, j));
			dyy = one_over_dy2*(u(i, j + 1) - 2*u(i, j) + u(i, j - 1));

			// v, vp derivatives
			dv    = (T)1/sqrt(2)*one_over_2dx*(u(i + 1, j + 1) - u(i - 1, j - 1));
			dvp   = (T)1/sqrt(2)*one_over_2dy*(u(i + 1, j - 1) - u(i - 1, j + 1));
			dvv   = (T)1/2*one_over_dx2*(u(i + 1, j + 1) + u(i - 1, j - 1) - 2*u(i, j));
			dvpvp = (T)1/2*one_over_dx2*(u(i + 1, j - 1) + u(i - 1, j + 1) - 2*u(i, j));

			// MA_1 and MA_2 - Need to add the constant term u_0
			if (is_regularized)
			{
				(MA_array.values[0])(i, j) = (T)0.5*(dxx + sqrt(POW2(dxx) + POW2(delta)))*(T)0.5*(dyy + sqrt(POW2(dyy) + POW2(delta))) + (T)0.5*(dxx - sqrt(POW2(dxx) + POW2(delta))) + (T)0.5*(dyy - sqrt(POW2(dyy) + POW2(delta))) - rho_x(i, j)/rho_y_1(i, j) - u(mid_i, mid_j);
				(MA_array.values[1])(i, j) = (T)0.5*(dvv + sqrt(POW2(dvv) + POW2(delta)))*(T)0.5*(dvpvp + sqrt(POW2(dvpvp) + POW2(delta))) + (T)0.5*(dvv - sqrt(POW2(dvv) + POW2(delta))) + (T)0.5*(dvpvp - sqrt(POW2(dvpvp) + POW2(delta))) - rho_x(i, j)/rho_y_2(i, j) - u(mid_i, mid_j);
			}
			else
			{
				(MA_array.values[0])(i, j) = max(dxx, delta)*max(dyy, delta) + min(dxx, delta) + min(dyy, delta) - rho_x(i, j)/rho_y_1(i, j) - u(mid_i, mid_j);
				(MA_array.values[1])(i, j) = max(dvv, delta)*max(dvpvp, delta) + min(dvv, delta) + min(dvpvp, delta) - rho_x(i, j)/rho_y_2(i, j) - u(mid_i, mid_j);
			}
						
			if (is_regularized)
			{
				MA(i, j) = (T)0.5*((MA_array.values[0])(i, j) + (MA_array.values[1])(i, j) + sqrt(POW2((MA_array.values[0])(i, j) - (MA_array.values[1])(i, j)) + POW2(delta)));
			}
			else
			{
				MA(i, j) = min((MA_array.values[0])(i, j), (MA_array.values[1])(i, j));
			}
		}
		END_GRID_ITERATION_2D;
	}

	void CalculateJacobian(CSR_MATRIX<T>& J, VECTOR_ND<T>& b_vector, const int& thread_id)
	{
		AssignStencil(density_x, density_y, density_y, thread_id);
		
		const int num_all_full_cells = AssignSequentialindicesToFullCells(bc, thread_id);
		const int nnz = CountNonZeroElements(bc, thread_id, MA_array, num_all_full_cells);
		
		HEAD_THREAD_WORK(J.Initialize(num_all_full_cells, nnz, multithreading));
		HEAD_THREAD_WORK(b_vector.Initialize(num_all_full_cells));
		HEAD_THREAD_WORK(multithreading->SplitDomainIndex1D(0, num_all_full_cells));

		BEGIN_HEAD_THREAD_WORK
		{
			J.start_ix[0] = 0;
			J.end_ix[0] = multithreading->sync_value_int[0] - 1;
			J.prev_row_array[0] = -1;
			J.values_ix_array[0] = 0;
			for (int thread_id = 1; thread_id < multithreading->num_threads; thread_id++)
			{
				J.start_ix[thread_id] = J.end_ix[thread_id - 1] + 1;
				J.end_ix[thread_id] = J.end_ix[thread_id - 1] + multithreading->sync_value_int[thread_id];
				J.prev_row_array[thread_id] = -1;
				J.values_ix_array[thread_id] = J.start_ix[thread_id];
			}
		}
		END_HEAD_THREAD_WORK

		// Speedup Variable
		T one_over_dx  = u.grid.one_over_dx, one_over_dy = u.grid.one_over_dy;
		T one_over_dx2 = u.grid.one_over_dx2, one_over_dy2 = u.grid.one_over_dy2;
		T one_over_2dx = u.grid.one_over_2dx, one_over_2dy = u.grid.one_over_2dy;
		
		T x_min = u.grid.x_min, x_max = u.grid.x_max, dx = u.grid.dx;
		T y_min = u.grid.y_min, y_max = u.grid.y_max, dy = u.grid.dy;

		int i_start = u.grid.i_start, i_end = u.grid.i_end;
		int j_start = u.grid.j_start, j_end = u.grid.j_end;
		
		// For fixed point
		int mid_i(0), mid_j(0);
		
		mid_i = (T)0.5*(i_start + i_end);
		mid_j = (T)0.5*(j_start + j_end);

		T H_star_n(0);

		ARRAY<VT> target_normal;
		target_normal.Initialize(n_y);

		for (int i = 0; i < n_y; i++)
		{
			target_normal[i].x = cos(2 * PI*i / n_y);
			target_normal[i].y = sin(2 * PI*i / n_y);
		}

		// Define H*(n)
		H_n.Initialize(n_y);

		if (test_number == 1)
		{
			for (int k = 0; k < target_normal.length; k++)
			{
				VT& tn(target_normal[k]);

				T x_coor_0 = x_min + i_start*dx, y_coor_0 = y_min + j_start*dy;

				VT y_start(x_coor_0, y_coor_0);

				H_n[k] = DotProduct(y_start, tn);

				// Left Boundary
				for (int j = j_start; j <= j_end; j++)
				{
					T x_coor = x_min + i_start*dx, y_coor = y_min + j*dy;

					VT y_0(x_coor, y_coor);

					H_n[k] = MAX(H_n[k], DotProduct(y_0, tn));
				}

				// Right Boundary
				for (int j = j_start; j <= j_end; j++)
				{
					T x_coor = x_min + i_end*dx, y_coor = y_min + j*dy;

					VT y_0(x_coor, y_coor);

					H_n[k] = MAX(H_n[k], DotProduct(y_0, tn));
				}

				// Lower Boundary
				for (int i = i_start; i <= i_end; i++)
				{
					T x_coor = x_min + i*dx, y_coor = y_min + j_start*dy;

					VT y_0(x_coor, y_coor);

					H_n[k] = MAX(H_n[k], DotProduct(y_0, tn));
				}

				// Upper Boundary
				for (int i = i_start; i <= i_end; i++)
				{
					T x_coor = x_min + i*dx, y_coor = y_min + j_end*dy;

					VT y_0(x_coor, y_coor);

					H_n[k] = MAX(H_n[k], DotProduct(y_0, tn));
				} 
			} 
		}

		if (test_number == 2)
		{
			/*int num_of_target_boundary(0);

			BEGIN_GRID_ITERATION_2D(u.partial_grids[thread_id])
			{
				T x_coor = u.grid.x_min + i*u.grid.dx, y_coor = u.grid.y_min + j*u.grid.dy;
				T x_coor_l = u.grid.x_min + (i - 1)*u.grid.dx, y_coor_d = u.grid.y_min + (j - 1)*u.grid.dy;
				T x_coor_r = u.grid.x_min + (i + 1)*u.grid.dx, y_coor_u = u.grid.y_min + (j + 1)*u.grid.dy;
				
				T domain_region = POW2(x_coor)/POW2(0.8) + POW2(y_coor)/POW2(0.4) - 1;
				T domain_region_r = POW2(x_coor_r)/POW2(0.8) + POW2(y_coor)/POW2(0.4) - 1;
				T domain_region_l = POW2(x_coor_l)/POW2(0.8) + POW2(y_coor)/POW2(0.4) - 1;
				T domain_region_u = POW2(x_coor)/POW2(0.8) + POW2(y_coor_u)/POW2(0.4) - 1;
				T domain_region_d = POW2(x_coor)/POW2(0.8) + POW2(y_coor_d)/POW2(0.4) - 1;

				if (domain_region < 0)
				{
					if (domain_region_r > 0)
					{
						num_of_target_boundary++;
					}
					else if (domain_region_l > 0)
					{
						num_of_target_boundary++;
					}
					else if (domain_region_d > 0)
					{
						num_of_target_boundary++;
					}
					else if (domain_region_u > 0)
					{
						num_of_target_boundary++;
					}
				}
				if (domain_region == 0)
				{
					num_of_target_boundary++;
				}
			}
			END_GRID_ITERATION_2D;
			
			ARRAY<VT> target_boundary;
			target_boundary.Initialize(num_of_target_boundary);
			
			int index_for_target_boundary(0);

			BEGIN_GRID_ITERATION_2D(u.partial_grids[thread_id])
			{
				T x_coor = u.grid.x_min + i*u.grid.dx, y_coor = u.grid.y_min + j*u.grid.dy;
				T x_coor_l = u.grid.x_min + (i - 1)*u.grid.dx, y_coor_d = u.grid.y_min + (j - 1)*u.grid.dy;
				T x_coor_r = u.grid.x_min + (i + 1)*u.grid.dx, y_coor_u = u.grid.y_min + (j + 1)*u.grid.dy;

				T domain_region = POW2(x_coor)/POW2(0.8) + POW2(y_coor)/POW2(0.4) - 1;
				T domain_region_r = POW2(x_coor_r)/POW2(0.8) + POW2(y_coor)/POW2(0.4) - 1;
				T domain_region_l = POW2(x_coor_l)/POW2(0.8) + POW2(y_coor)/POW2(0.4) - 1;
				T domain_region_u = POW2(x_coor)/POW2(0.8) + POW2(y_coor_u)/POW2(0.4) - 1;
				T domain_region_d = POW2(x_coor)/POW2(0.8) + POW2(y_coor_d)/POW2(0.4) - 1;

				MATRIX_2X2 M_x(0.8, 0, 0, 0.4), M_y(0.6, 0.2, 0.2, 0.8), J(0, 1, -1, 0);

				MATRIX_2X2 nu, deno;

				nu = M_x.Inversed()*M_y.Inversed()*J;
				deno = M_x.Inversed()*M_y.Inversed();
			
				T atheta = nu.Trace()/deno.Trace();

				MATRIX_2X2 R_theta(cos(atan(atheta)), sin(atan(atheta)), -sin(atan(atheta)), cos(atan(atheta)));

				MATRIX_2X2 Ellipse_Mapping;

				Ellipse_Mapping = M_y*R_theta*M_x.Inversed();

				T x_vec_c = Ellipse_Mapping.x[0]*x_coor + Ellipse_Mapping.x[2]*y_coor, y_vec_c = Ellipse_Mapping.x[1]*x_coor + Ellipse_Mapping.x[3]*y_coor;

				if (domain_region < 0)
				{
					if (domain_region_r > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_vec_c, y_vec_c);
						index_for_target_boundary++;
					}
					else if (domain_region_l > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_vec_c, y_vec_c);
						index_for_target_boundary++;
					}
					else if (domain_region_d > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_vec_c, y_vec_c);
						index_for_target_boundary++;
					}
					else if (domain_region_u > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_vec_c, y_vec_c);
						index_for_target_boundary++;
					}
				}
				if (domain_region == 0)
				{
					target_boundary[index_for_target_boundary] = VT(x_vec_c, y_vec_c);
					index_for_target_boundary++;
				}
			}
			END_GRID_ITERATION_2D;

			for (int k = 0; k < n_y; k++)
			{
				VT& tn(target_normal[k]);

				H_n[k] = DotProduct(target_boundary[0], tn);
				
				for (int i = 0; i < num_of_target_boundary; i++)
				{
					H_n[k] = MAX(H_n[k], DotProduct(target_boundary[i], tn));
				}
			}*/
			
			// Boundary Vector define
			MATRIX_2X2 M_x(0.8, 0, 0, 0.4), M_y(0.6, 0.2, 0.2, 0.8), J(0, 1, -1, 0);

			MATRIX_2X2 nu, deno;

			nu = M_x.Inversed()*M_y.Inversed()*J;
			deno = M_x.Inversed()*M_y.Inversed();
			
			T atheta = nu.Trace()/deno.Trace();

			MATRIX_2X2 R_theta(cos(atan(atheta)), sin(atan(atheta)), -sin(atan(atheta)), cos(atan(atheta)));

			MATRIX_2X2 Ellipse_Mapping;

			Ellipse_Mapping = M_y*R_theta*M_x.Inversed();

			ARRAY<VT> target_boundary;
			target_boundary.Initialize(N_point);
				
			for (int i = 0; i < N_point; i++)
			{
				T angle = 2*i*PI/N_point;

				target_boundary[i].x = 0.8*Ellipse_Mapping.x[0]*cos(angle) + 0.4*Ellipse_Mapping.x[2]*sin(angle);
				target_boundary[i].y = 0.8*Ellipse_Mapping.x[1]*cos(angle) + 0.4*Ellipse_Mapping.x[3]*sin(angle);

				/*target_boundary[i].x = 0.6*cos(angle) - 0.2*sin(angle);
				target_boundary[i].y = 0.2*cos(angle) - 0.8*sin(angle);*/
			}
				
			for (int k = 0; k < n_y; k++)
			{
				VT& tn(target_normal[k]);

				H_n[k] = DotProduct(target_boundary[0], tn);
				
				for (int i = 0; i < N_point; i++)
				{
					H_n[k] = MAX(H_n[k], DotProduct(target_boundary[i], tn));
				}
			}
		}

		if (test_number == 3)
		{
			int num_of_target_boundary(0);

			BEGIN_GRID_ITERATION_2D(u.partial_grids[thread_id])
			{
				T x_coor = u.grid.x_min + i*u.grid.dx, y_coor = u.grid.y_min + j*u.grid.dy;
				T x_coor_l = u.grid.x_min + (i - 1)*u.grid.dx, y_coor_d = u.grid.y_min + (j - 1)*u.grid.dy;
				T x_coor_r = u.grid.x_min + (i + 1)*u.grid.dx, y_coor_u = u.grid.y_min + (j + 1)*u.grid.dy;

				T domain_region_1 = POW2(x_coor + 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_1_r = POW2(x_coor_r + 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_1_l = POW2(x_coor_l + 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_1_u = POW2(x_coor + 0.1) + POW2(y_coor_u) - POW2(0.85);
				T domain_region_1_d = POW2(x_coor + 0.1) + POW2(y_coor_d) - POW2(0.85);

				T domain_region_2 = POW2(x_coor - 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_2_r = POW2(x_coor_r - 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_2_l = POW2(x_coor_l - 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_2_u = POW2(x_coor - 0.1) + POW2(y_coor_u) - POW2(0.85);
				T domain_region_2_d = POW2(x_coor - 0.1) + POW2(y_coor_d) - POW2(0.85);

				if (domain_region_1 < 0 && x_coor < -0.1)
				{
					if (domain_region_1_r > 0)
					{
						num_of_target_boundary++;
					}
					else if (domain_region_1_l > 0)
					{
						num_of_target_boundary++;
					}
					else if (domain_region_1_d > 0)
					{
						num_of_target_boundary++;
					}
					else if (domain_region_1_u > 0)
					{
						num_of_target_boundary++;
					}
				}
				if (domain_region_2 < 0 && x_coor > 0.1)
				{
					if (domain_region_2_r > 0)
					{
						num_of_target_boundary++;
					}
					else if (domain_region_2_l > 0)
					{
						num_of_target_boundary++;
					}
					else if (domain_region_2_d > 0)
					{
						num_of_target_boundary++;
					}
					else if (domain_region_2_u > 0)
					{
						num_of_target_boundary++;
					}
				}
			}
			END_GRID_ITERATION_2D;
			
			ARRAY<VT> target_boundary;
			target_boundary.Initialize(num_of_target_boundary);
			
			int index_for_target_boundary(0);

			BEGIN_GRID_ITERATION_2D(u.partial_grids[thread_id])
			{
				T x_coor = u.grid.x_min + i*u.grid.dx, y_coor = u.grid.y_min + j*u.grid.dy;
				T x_coor_l = u.grid.x_min + (i - 1)*u.grid.dx, y_coor_d = u.grid.y_min + (j - 1)*u.grid.dy;
				T x_coor_r = u.grid.x_min + (i + 1)*u.grid.dx, y_coor_u = u.grid.y_min + (j + 1)*u.grid.dy;

				T domain_region_1 = POW2(x_coor + 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_1_r = POW2(x_coor_r + 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_1_l = POW2(x_coor_l + 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_1_u = POW2(x_coor + 0.1) + POW2(y_coor_u) - POW2(0.85);
				T domain_region_1_d = POW2(x_coor + 0.1) + POW2(y_coor_d) - POW2(0.85);

				T domain_region_2 = POW2(x_coor - 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_2_r = POW2(x_coor_r - 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_2_l = POW2(x_coor_l - 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_2_u = POW2(x_coor - 0.1) + POW2(y_coor_u) - POW2(0.85);
				T domain_region_2_d = POW2(x_coor - 0.1) + POW2(y_coor_d) - POW2(0.85);

				if (domain_region_1 < 0 && x_coor < -0.1)
				{
					if (domain_region_1_r > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_coor + 0.1, y_coor);
						index_for_target_boundary++;
					}
					else if (domain_region_1_l > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_coor + 0.1, y_coor);
						index_for_target_boundary++;
					}
					else if (domain_region_1_d > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_coor + 0.1, y_coor);
						index_for_target_boundary++;
					}
					else if (domain_region_1_u > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_coor + 0.1, y_coor);
						index_for_target_boundary++;
					}
				}
				if (domain_region_2 < 0 && x_coor > 0.1)
				{
					if (domain_region_2_r > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_coor - 0.1, y_coor);
						index_for_target_boundary++;
					}
					else if (domain_region_2_l > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_coor - 0.1, y_coor);
						index_for_target_boundary++;
					}
					else if (domain_region_2_d > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_coor - 0.1, y_coor);
						index_for_target_boundary++;
					}
					else if (domain_region_2_u > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_coor - 0.1, y_coor);
						index_for_target_boundary++;
					}
				}
			}
			END_GRID_ITERATION_2D;

			for (int k = 0; k < n_y; k++)
			{
				VT& tn(target_normal[k]);

				H_n[k] = DotProduct(target_boundary[0], tn);
				
				for (int i = 0; i < num_of_target_boundary; i++)
				{
					H_n[k] = MAX(H_n[k], DotProduct(target_boundary[i], tn));
				}
			}

			// Boundary Vector define
			/*ARRAY<VT> target_boundary;
			target_boundary.Initialize(N_point);
				
			for (int i = 0; i < N_point; i++)
			{
				T angle = 2*i*PI/N_point;

				target_boundary[i].x = 0.85*cos(angle);
				target_boundary[i].y = 0.85*sin(angle);
			}
				
			for (int k = 0; k < n_y; k++)
			{
				VT& tn(target_normal[k]);

				H_n[k] = DotProduct(target_boundary[0], tn);
				
				for (int i = 0; i < N_point; i++)
				{
					H_n[k] = MAX(H_n[k], DotProduct(target_boundary[i], tn));
				}
			}*/
		}

		if (test_number == 8)
		{
			// Boundary Vector define
			T a = cos(PI/8), radius = sqrt((T)1/POW2(a) - 1);
			
			ARRAY<VT> target_boundary;
			target_boundary.Initialize(N_point);
				
			for (int i = 0; i < N_point; i++)
			{
				target_boundary[i].x = radius*cos(2*i*PI/N_point);
				target_boundary[i].y = radius*sin(2*i*PI/N_point) + (T)1/a;
			}
				
			for (int k = 0; k < n_y; k++)
			{
				VT& tn(target_normal[k]);

				H_n[k] = DotProduct(target_boundary[0], tn);
				
				for (int i = 0; i < N_point; i++)
				{
					H_n[k] = MAX(H_n[k], DotProduct(target_boundary[i], tn));
				}
			}
		}
				
		BEGIN_GRID_ITERATION_2D(u.partial_grids[thread_id])
		{
			T dxx, dyy, dvv, dvpvp;

			dxx = one_over_dx2*(u(i + 1, j) - 2*u(i, j) + u(i - 1, j));
			dyy = one_over_dy2*(u(i, j + 1) - 2*u(i, j) + u(i, j - 1));
			dvv   = (T)1/2*one_over_dx2*(u(i + 1, j + 1) + u(i - 1, j - 1) - 2*u(i, j));
			dvpvp = (T)1/2*one_over_dx2*(u(i + 1, j - 1) + u(i - 1, j + 1) - 2*u(i, j));

			if (bc(i, j) < 0)
			{
				continue;
			}
						
			// Need to add the derivative of source term
			if (bc(i - 1, j) < 0)
			{
				if (bc(i, j - 1) < 0)
				{
					// Define b_vector
					T dxp = one_over_dx*(u(i + 1, j) - u(i, j)), dyp = one_over_dy*(u(i, j + 1) - u(i, j));

					T b_vector_value = 0;
					int max_left_lower_corner(0);
					int sub_k(0);

					for (int k = 0; k < target_normal.length; k++)
					{
						VT& tn(target_normal[k]);

						T angle = 2*PI*k/n_y;
												
						if (angle >= PI && angle <= 3*PI/2)
						{
							if (sub_k == 0)
							{
								b_vector_value = tn.x*dxp + tn.y*dyp - H_n[k];	
							}

							T compared_one = tn.x*dxp + tn.y*dyp - H_n[k];	
							
							if (b_vector_value <= compared_one)
							{
								b_vector_value = compared_one;
								max_left_lower_corner = k;
							}
							
							sub_k++;
						}
						else
						{
							continue;
						}
					}

					b_vector[bc(i, j)] = b_vector_value;

					T u_r, u_u, u_c;
										
					u_r = target_normal[max_left_lower_corner].x*one_over_dx;
					u_u = target_normal[max_left_lower_corner].y*one_over_dy;
					u_c = -target_normal[max_left_lower_corner].x*one_over_dx - target_normal[max_left_lower_corner].y*one_over_dy;

					if (bc(i + 1, j) > -1)
					{
						J.AssignValue(bc(i, j), bc(i + 1, j), u_r, thread_id);
					}
					if (bc(i, j + 1) > -1)
					{
						J.AssignValue(bc(i, j), bc(i, j + 1), u_u, thread_id);
					}	
			
					J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
				}
				else if (bc(i, j + 1) < 0)
				{
					// Define b_vector
					T dxp = one_over_dx*(u(i + 1, j) - u(i, j)), dym = one_over_dy*(u(i, j) - u(i, j - 1));

					T b_vector_value = 0;
					int max_left_upper_corner(0);
					int sub_k(0);

					for (int k = 0; k < target_normal.length; k++)
					{
						VT& tn(target_normal[k]);

						T angle = 2*PI*k/n_y;
												
						if (angle >= PI/2 && angle <= PI)
						{
							if (sub_k == 0)
							{
								b_vector_value = tn.x*dxp + tn.y*dym - H_n[k];	
							}

							T compared_one = tn.x*dxp + tn.y*dym - H_n[k];	
							
							if (b_vector_value <= compared_one)
							{
								b_vector_value = compared_one;
								max_left_upper_corner = k;
							}

							sub_k++;
						}
						else
						{
							continue;
						}
					}

					b_vector[bc(i, j)] = b_vector_value;

					T u_r, u_d, u_c;

					u_r = target_normal[max_left_upper_corner].x*one_over_dx;
					u_d = -target_normal[max_left_upper_corner].y*one_over_dy;
					u_c = -target_normal[max_left_upper_corner].x*one_over_dx + target_normal[max_left_upper_corner].y*one_over_dy;

					if (bc(i + 1, j) > -1)
					{
						J.AssignValue(bc(i, j), bc(i + 1, j), u_r, thread_id);
					}
					if (bc(i, j - 1) > -1)
					{
						J.AssignValue(bc(i, j), bc(i, j - 1), u_d, thread_id);
					}	

					J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
				}
				else
				{
					// Define b_vector
					T dxp = one_over_dx*(u(i + 1, j) - u(i, j)), dym = one_over_dy*(u(i, j) - u(i, j - 1)), dyp = one_over_dy*(u(i, j + 1) - u(i, j));

					T b_vector_value(0);
					int sub_k(0);
					int max_left(0);

					for (int k = 0; k < target_normal.length; k++)
					{
						VT n_x(-1, 0);

						VT& tn(target_normal[k]);

						T dot = DotProduct(n_x, tn);

						if (dot <= 0)
						{
							continue;
						}
						else
						{
							if (sub_k == 0)
							{
								b_vector_value = tn.x*dxp + MAX(tn.y, 0)*dym + MIN(tn.y, 0)*dyp - H_n[k];
							}

							T compared_one = tn.x*dxp + MAX(tn.y, 0)*dym + MIN(tn.y, 0)*dyp - H_n[k];
						
							if (b_vector_value <= compared_one)
							{
								b_vector_value = compared_one;
								max_left = k;
							}
						
							sub_k++;
						}
					}

					b_vector[bc(i, j)] = b_vector_value;

					if (target_normal[max_left].y > 0)
					{
						T u_r, u_d, u_c;

						u_r = target_normal[max_left].x*one_over_dx;
						u_d = -target_normal[max_left].y*one_over_dy;
						u_c = -target_normal[max_left].x*one_over_dx + target_normal[max_left].y*one_over_dy;

						if (bc(i + 1, j) > -1)
						{
							J.AssignValue(bc(i, j), bc(i + 1, j), u_r, thread_id);
						}
						if (bc(i, j - 1) > -1)
						{
							J.AssignValue(bc(i, j), bc(i, j - 1), u_d, thread_id);
						}
						
						J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
					}
					else if (target_normal[max_left].y < 0)
					{
						T u_r, u_u, u_c;

						u_r = target_normal[max_left].x*one_over_dx;
						u_u = target_normal[max_left].y*one_over_dy;
						u_c = -target_normal[max_left].x*one_over_dx - target_normal[max_left].y*one_over_dy;

						if (bc(i + 1, j) > -1)
						{
							J.AssignValue(bc(i, j), bc(i + 1, j), u_r, thread_id);
						}
						if (bc(i, j + 1) > -1)
						{
							J.AssignValue(bc(i, j), bc(i, j + 1), u_u, thread_id);
						}

						J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
					}
					else
					{
						T u_r, u_c;

						u_r = target_normal[max_left].x*one_over_dx;
						u_c = -target_normal[max_left].x*one_over_dx;

						if (bc(i + 1, j) > -1)
						{
							J.AssignValue(bc(i, j), bc(i + 1, j), u_r, thread_id);
						}
						
						J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
					}
				}
			}
			else if (bc(i + 1, j) < 0)
			{
				if (bc(i, j - 1) < 0)
				{
					// Define b_vector
					T dxm = one_over_dx*(u(i, j) - u(i - 1, j)), dyp = one_over_dy*(u(i, j + 1) - u(i, j));

					T b_vector_value = 0;
					int max_right_lower_corner(0);
					int sub_k(0);

					for (int k = 0; k < target_normal.length; k++)
					{
						VT& tn(target_normal[k]);

						T angle = 2*PI*k/n_y;
												
						if (angle >= 3*PI/2 && angle <= 2*PI)
						{
							if (sub_k == 0)
							{
								b_vector_value = tn.x*dxm + tn.y*dyp - H_n[k];	
							}

							T compared_one = tn.x*dxm + tn.y*dyp - H_n[k];	
							
							if (b_vector_value <= compared_one)
							{
								b_vector_value = compared_one;
								max_right_lower_corner = k;
							}
							sub_k++;
						}
						else
						{
							continue;
						}
					}

					b_vector[bc(i, j)] = b_vector_value;

					T u_l, u_u, u_c;

					u_l = -target_normal[max_right_lower_corner].x*one_over_dx;
					u_u = target_normal[max_right_lower_corner].y*one_over_dy;
					u_c = target_normal[max_right_lower_corner].x*one_over_dx - target_normal[max_right_lower_corner].y*one_over_dy;

					if (bc(i - 1, j) > -1)
					{
						J.AssignValue(bc(i, j), bc(i - 1, j), u_l, thread_id);
					}
					if (bc(i, j + 1) > -1)
					{
						J.AssignValue(bc(i, j), bc(i, j + 1), u_u, thread_id);
					}	
				
					J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
				}
				else if (bc(i, j + 1) < 0)
				{
					// Define b_vector
					T dxm = one_over_dx*(u(i, j) - u(i - 1, j)), dym = one_over_dy*(u(i, j) - u(i, j - 1));

					T b_vector_value = 0;
					int max_right_upper_corner(0);
					int sub_k(0);

					for (int k = 0; k < target_normal.length; k++)
					{
						VT& tn(target_normal[k]);

						T angle = 2*PI*k/n_y;
												
						if (angle >= 0 && angle <= PI/2)
						{
							if (sub_k ==0)
							{
								b_vector_value = tn.x*dxm + tn.y*dym - H_n[k];	
							}
							
							T compared_one = tn.x*dxm + tn.y*dym - H_n[k];	
							
							if (b_vector_value <= compared_one)
							{
								b_vector_value = compared_one;
								max_right_upper_corner = k;
							}
							sub_k++;
						}
						else
						{
							continue;
						}
					}

					b_vector[bc(i, j)] = b_vector_value;

					T u_l, u_d, u_c;

					u_l = -target_normal[max_right_upper_corner].x*one_over_dx;
					u_d = -target_normal[max_right_upper_corner].y*one_over_dy;
					u_c = target_normal[max_right_upper_corner].x*one_over_dx + target_normal[max_right_upper_corner].y*one_over_dy;

					if (bc(i - 1, j) > -1)
					{
						J.AssignValue(bc(i, j), bc(i - 1, j), u_l, thread_id);
					}
					if (bc(i, j - 1) > -1)
					{
						J.AssignValue(bc(i, j), bc(i, j - 1), u_d, thread_id);
					}	
				
					J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
				}
				else
				{
					// Define b_vector
					T dxm = one_over_dx*(u(i, j) - u(i - 1, j)), dym = one_over_dy*(u(i, j) - u(i, j - 1)), dyp = one_over_dy*(u(i, j + 1) - u(i, j));

					T b_vector_value(0);
					int max_right(0);
					int sub_k(0);

					for (int k = 0; k < target_normal.length; k++)
					{
						VT n_x(1, 0);
						VT& tn(target_normal[k]);

						T dot = DotProduct(n_x, tn);

						if (dot <= 0)
						{
							continue;
						}

						if (sub_k == 0)
						{
							b_vector_value = tn.x*dxm + MAX(tn.y, 0)*dym + MIN(tn.y, 0)*dyp - H_n[k];
						}

						T compared_one = tn.x*dxm + MAX(tn.y, 0)*dym + MIN(tn.y, 0)*dyp - H_n[k];

						if (b_vector_value <= compared_one)
						{
							b_vector_value = compared_one;
							max_right = k;
						}
						sub_k++;
					}

					b_vector[bc(i, j)] = b_vector_value;

					if (target_normal[max_right].y > 0)
					{
						T u_l, u_d, u_c;

						u_l = -target_normal[max_right].x*one_over_dx;
						u_d = -target_normal[max_right].y*one_over_dy;
						u_c = target_normal[max_right].x*one_over_dx + target_normal[max_right].y*one_over_dy;

						if (bc(i - 1, j) > -1)
						{
							J.AssignValue(bc(i, j), bc(i - 1, j), u_l, thread_id);
						}
						if (bc(i, j - 1) > -1)
						{
							J.AssignValue(bc(i, j), bc(i, j - 1), u_d, thread_id);
						}
						
						J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
					}
					else if (target_normal[max_right].y < 0)
					{
						T u_l, u_u, u_c;

						u_l = -target_normal[max_right].x*one_over_dx;
						u_u = target_normal[max_right].y*one_over_dy;
						u_c = target_normal[max_right].x*one_over_dx - target_normal[max_right].y*one_over_dy;

						if (bc(i - 1, j) > -1)
						{
							J.AssignValue(bc(i, j), bc(i - 1, j), u_l, thread_id);
						}
						if (bc(i, j + 1) > -1)
						{
							J.AssignValue(bc(i, j), bc(i, j + 1), u_u, thread_id);
						}

						J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
					}
					else
					{
						T u_l, u_c;

						u_l = -target_normal[max_right].x*one_over_dx;
						u_c = target_normal[max_right].x*one_over_dx;

						if (bc(i - 1, j) > -1)
						{
							J.AssignValue(bc(i, j), bc(i - 1, j), u_l, thread_id);
						}
						
						J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
					}
				}
			}
			else if (bc(i, j - 1) < 0)
			{
				// Define b_vector
				T dxm = one_over_dx*(u(i, j) - u(i - 1, j)), dxp = one_over_dx*(u(i + 1, j) - u(i, j)), dyp = one_over_dy*(u(i, j + 1) - u(i, j));

				T b_vector_value(0);
				int sub_k(0);

				for (int k = 0; k < target_normal.length; k++)
				{
					VT n_y(0, -1);
					VT& tn(target_normal[k]);

					T dot = DotProduct(n_y, tn);

					if (dot <= 0)
					{
						continue;
					}

					if (sub_k == 0)
					{
						b_vector_value = MAX(tn.x, 0)*dxm + MIN(tn.x, 0)*dxp + tn.y*dyp - H_n[k];
					}

					T compared_one = MAX(tn.x, 0)*dxm + MIN(tn.x, 0)*dxp + tn.y*dyp - H_n[k];

					if (b_vector_value <= compared_one)
					{
						b_vector_value = compared_one;
						max_lower = k;
					}
					sub_k++;
				}

				b_vector[bc(i, j)] = b_vector_value;

				if (target_normal[max_lower].x > 0)
				{
					T u_l, u_u, u_c;

					u_l = -target_normal[max_lower].x*one_over_dx;
					u_u = target_normal[max_lower].y*one_over_dy;
					u_c = target_normal[max_lower].x*one_over_dx - target_normal[max_lower].y*one_over_dy;

					if (bc(i - 1, j) > -1)
					{
						J.AssignValue(bc(i, j), bc(i - 1, j), u_l, thread_id);
					}
					if (bc(i, j + 1) > -1)
					{
						J.AssignValue(bc(i, j), bc(i, j + 1), u_u, thread_id);
					}
						
					J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
				}
				else if (target_normal[max_lower].x < 0)
				{
					T u_r, u_u, u_c;

					u_r = target_normal[max_lower].x*one_over_dx;
					u_u = target_normal[max_lower].y*one_over_dy;
					u_c = -target_normal[max_lower].x*one_over_dx - target_normal[max_lower].y*one_over_dy;

					if (bc(i + 1, j) > -1)
					{
						J.AssignValue(bc(i, j), bc(i + 1, j), u_r, thread_id);
					}
					if (bc(i, j + 1) > -1)
					{
						J.AssignValue(bc(i, j), bc(i, j + 1), u_u, thread_id);
					}

					J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
				}
				else
				{
					T u_u, u_c;

					u_u = target_normal[max_lower].y*one_over_dy;
					u_c = -target_normal[max_lower].y*one_over_dy;

					if (bc(i, j + 1) > -1)
					{
						J.AssignValue(bc(i, j), bc(i, j + 1), u_u, thread_id);
					}
						
					J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
				}
			}
			else if (bc(i, j + 1) < 0)
			{
				// Define b_vector
				T dxm = one_over_dx*(u(i, j) - u(i - 1, j)), dxp = one_over_dx*(u(i + 1, j) - u(i, j)), dym = one_over_dy*(u(i, j) - u(i, j - 1));

				T b_vector_value = 0;
				int sub_k(0);

				for (int k = 0; k < target_normal.length; k++)
				{
					VT n_y(0, 1);
					VT& tn(target_normal[k]);

					T dot = DotProduct(n_y, tn);

					if (dot <= 0)
					{
						continue;
					}

					if (sub_k == 0)
					{
						b_vector_value = MAX(tn.x, 0)*dxm + MIN(tn.x, 0)*dxp + tn.y*dym - H_n[k];
					}

					T compared_one = MAX(tn.x, 0)*dxm + MIN(tn.x, 0)*dxp + tn.y*dym - H_n[k];

					if (b_vector_value <= compared_one)
					{
						b_vector_value = compared_one;
						max_upper = k;
					}
					sub_k++;
				}

				b_vector[bc(i, j)] = b_vector_value;

				if (target_normal[max_upper].x > 0)
				{
					T u_l, u_d, u_c;

					u_l = -target_normal[max_upper].x*one_over_dx;
					u_d = -target_normal[max_upper].y*one_over_dy;
					u_c = target_normal[max_upper].x*one_over_dx + target_normal[max_upper].y*one_over_dy;

					if (bc(i - 1, j) > -1)
					{
						J.AssignValue(bc(i, j), bc(i - 1, j), u_l, thread_id);
					}
					if (bc(i, j - 1) > -1)
					{
						J.AssignValue(bc(i, j), bc(i, j - 1), u_d, thread_id);
					}
						
					J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
				}
				else if (target_normal[max_upper].x < 0)
				{
					T u_r, u_d, u_c;

					u_r = target_normal[max_upper].x*one_over_dx;
					u_d = -target_normal[max_upper].y*one_over_dy;
					u_c = -target_normal[max_upper].x*one_over_dx + target_normal[max_upper].y*one_over_dy;

					if (bc(i + 1, j) > -1)
					{
						J.AssignValue(bc(i, j), bc(i + 1, j), u_r, thread_id);
					}
					if (bc(i, j - 1) > -1)
					{
						J.AssignValue(bc(i, j), bc(i, j - 1), u_d, thread_id);
					}

					J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
				}
				else
				{
					T u_d, u_c;

					u_d = -target_normal[max_upper].y*one_over_dy;
					u_c = target_normal[max_upper].y*one_over_dy;

					if (bc(i, j - 1) > -1)
					{
						J.AssignValue(bc(i, j), bc(i, j - 1), u_d, thread_id);
					}
						
					J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
				}
			}
			else
			{
				if (is_regularized)
				{
					b_vector[bc(i, j)] = MA(i, j);

					T max_dxx_r, max_dyy_r, max_dvv_r, max_dvpvp_r, min_dxx_r, min_dyy_r, min_dvv_r, min_dvpvp_r;

					max_dxx_r = (T)0.5*(dxx + sqrt(POW2(dxx) + POW2(delta)));
					max_dyy_r = (T)0.5*(dyy + sqrt(POW2(dyy) + POW2(delta)));
					min_dxx_r = (T)0.5*(dxx - sqrt(POW2(dxx) + POW2(delta)));
					min_dyy_r = (T)0.5*(dyy - sqrt(POW2(dyy) + POW2(delta)));

					max_dvv_r   = (T)0.5*(dvv + sqrt(POW2(dvv) + POW2(delta)));
					max_dvpvp_r = (T)0.5*(dvpvp + sqrt(POW2(dvpvp) + POW2(delta)));
					min_dvv_r   = (T)0.5*(dvv - sqrt(POW2(dvv) + POW2(delta)));
					min_dvpvp_r = (T)0.5*(dvpvp - sqrt(POW2(dvpvp) + POW2(delta)));

					T u_l, u_r, u_c, u_d, u_u, u_dl, u_dr, u_ul, u_ur;
										
					// Center Stencil
					T u_c_1, u_c_2, u_c_3;

					u_c_1 = 0;
					u_c_1 += (T)0.5*(-2*one_over_dx2 + dxx*(-2*one_over_dx2)/sqrt(POW2(dxx) + POW2(delta)))*max_dyy_r;
					u_c_1 += (T)0.5*(-2*one_over_dy2 + dyy*(-2*one_over_dy2)/sqrt(POW2(dyy) + POW2(delta)))*max_dxx_r;
					u_c_1 += (T)0.5*(-2*one_over_dx2 - dxx*(-2*one_over_dx2)/sqrt(POW2(dxx) + POW2(delta)));
					u_c_1 += (T)0.5*(-2*one_over_dy2 - dyy*(-2*one_over_dy2)/sqrt(POW2(dyy) + POW2(delta)));
					
					u_c_2 = 0;
					u_c_2 += (T)0.5*(-2*(T)0.5*one_over_dx2 + dvv*(-2*(T)0.5*one_over_dx2)/sqrt(POW2(dvv) + POW2(delta)))*max_dvpvp_r;
					u_c_2 += (T)0.5*(-2*(T)0.5*one_over_dy2 + dvpvp*(-2*(T)0.5*one_over_dy2)/sqrt(POW2(dvpvp) + POW2(delta)))*max_dvv_r;
					u_c_2 += (T)0.5*(-2*(T)0.5*one_over_dx2 - dvv*(-2*(T)0.5*one_over_dx2)/sqrt(POW2(dvv) + POW2(delta)));
					u_c_2 += (T)0.5*(-2*(T)0.5*one_over_dy2 - dvpvp*(-2*(T)0.5*one_over_dy2)/sqrt(POW2(dvpvp) + POW2(delta)));

					u_c_3 = 0;
					u_c_3 = (MA_array[0].array_for_this(i, j) - MA_array[1].array_for_this(i, j))*(u_c_1 - u_c_2)/sqrt(POW2(MA_array[0].array_for_this(i, j) - MA_array[1].array_for_this(i, j)) + POW2(delta));
										
					u_c = 0;
					
					u_c = u_c_1 + u_c_2 - u_c_3; 
					
					// Left Stencil
					T u_l_1, u_l_2, u_l_3;

					u_l_1 = 0;
					u_l_1 += (T)0.5*(one_over_dx2 + dxx*one_over_dx2/sqrt(POW2(dxx) + POW2(delta)))*max_dyy_r;
					u_l_1 += (T)0.5*(one_over_dx2 - dxx*one_over_dx2/sqrt(POW2(dxx) + POW2(delta)));
										
					u_l_2 = 0;
					
					u_l_3 = 0;
					u_l_3 = (MA_array[0].array_for_this(i, j) - MA_array[1].array_for_this(i, j))*(u_l_1 - u_l_2)/sqrt(POW2(MA_array[0].array_for_this(i, j) - MA_array[1].array_for_this(i, j)) + POW2(delta));
					
					u_l = 0;
					u_l = u_l_1 + u_l_2 - u_l_3; 
					
					
					// Right Stencil
					T u_r_1, u_r_2, u_r_3;

					u_r_1 = 0;
					u_r_1 += (T)0.5*(one_over_dx2 + dxx*one_over_dx2/sqrt(POW2(dxx) + POW2(delta)))*max_dyy_r;
					u_r_1 += (T)0.5*(one_over_dx2 - dxx*one_over_dx2/sqrt(POW2(dxx) + POW2(delta)));
									
					u_r_2 = 0;
					
					u_r_3 = 0;
					u_r_3 = (MA_array[0].array_for_this(i, j) - MA_array[1].array_for_this(i, j))*(u_r_1 - u_r_2)/sqrt(POW2(MA_array[0].array_for_this(i, j) - MA_array[1].array_for_this(i, j)) + POW2(delta));

					u_r = 0;
					u_r = u_r_1 + u_r_2 - u_r_3; 

					// Upper Stencil
					T u_u_1, u_u_2, u_u_3;

					u_u_1 = 0;
					u_u_1 += (T)0.5*(one_over_dy2 + dyy*one_over_dy2/sqrt(POW2(dyy) + POW2(delta)))*max_dxx_r;
					u_u_1 += (T)0.5*(one_over_dy2 - dyy*one_over_dy2/sqrt(POW2(dyy) + POW2(delta)));
					
					u_u_2 = 0;
					
					u_u_3 = 0;
					u_u_3 = (MA_array[0].array_for_this(i, j) - MA_array[1].array_for_this(i, j))*(u_u_1 - u_u_2)/sqrt(POW2(MA_array[0].array_for_this(i, j) - MA_array[1].array_for_this(i, j)) + POW2(delta));
					
					u_u = 0;
					u_u = u_u_1 + u_u_2 - u_u_3; 

					// Lower Stencil
					T u_d_1, u_d_2, u_d_3;

					u_d_1 = 0;
					u_d_1 += (T)0.5*(one_over_dy2 + dyy*one_over_dy2/sqrt(POW2(dyy) + POW2(delta)))*max_dxx_r;
					u_d_1 += (T)0.5*(one_over_dy2 - dyy*one_over_dy2/sqrt(POW2(dyy) + POW2(delta)));
					
					u_d_2 = 0;
					
					u_d_3 = 0;
					u_d_3 = (MA_array[0].array_for_this(i, j) - MA_array[1].array_for_this(i, j))*(u_d_1 - u_d_2)/sqrt(POW2(MA_array[0].array_for_this(i, j) - MA_array[1].array_for_this(i, j)) + POW2(delta));
					
					u_d = 0;
					u_d = u_d_1 + u_d_2 - u_d_3; 

					// Upper-Right Stencil
					T u_ur_1, u_ur_2, u_ur_3;

					u_ur_1 = 0;
										
					u_ur_2 = 0;
					u_ur_2 += (T)0.5*((T)0.5*one_over_dx2 + dvv*((T)0.5*one_over_dx2)/sqrt(POW2(dvv) + POW2(delta)))*max_dvpvp_r;
					u_ur_2 += (T)0.5*((T)0.5*one_over_dx2 - dvv*((T)0.5*one_over_dx2)/sqrt(POW2(dvv) + POW2(delta)));

					u_ur_3 = 0;
					u_ur_3 = (MA_array[0].array_for_this(i, j) - MA_array[1].array_for_this(i, j))*(u_ur_1 - u_ur_2)/sqrt(POW2(MA_array[0].array_for_this(i, j) - MA_array[1].array_for_this(i, j)) + POW2(delta));

					u_ur = 0;
					u_ur = u_ur_1 + u_ur_2 - u_ur_3; 

					// Down-Left Stencil
					T u_dl_1, u_dl_2, u_dl_3;

					u_dl_1 = 0;
										
					u_dl_2 = 0;
					u_dl_2 += (T)0.5*((T)0.5*one_over_dx2 + dvv*((T)0.5*one_over_dx2)/sqrt(POW2(dvv) + POW2(delta)))*max_dvpvp_r;
					u_dl_2 += (T)0.5*((T)0.5*one_over_dx2 - dvv*((T)0.5*one_over_dx2)/sqrt(POW2(dvv) + POW2(delta)));

					u_dl_3 = 0;
					u_dl_3 = (MA_array[0].array_for_this(i, j) - MA_array[1].array_for_this(i, j))*(u_dl_1 - u_dl_2)/sqrt(POW2(MA_array[0].array_for_this(i, j) - MA_array[1].array_for_this(i, j)) + POW2(delta));

					u_dl = 0;
					u_dl = u_dl_1 + u_dl_2 - u_dl_3; 
					
					// Upper-Left Stencil
					T u_ul_1, u_ul_2, u_ul_3;

					u_ul_1 = 0;
										
					u_ul_2 = 0;
					u_ul_2 += (T)0.5*((T)0.5*one_over_dy2 + dvpvp*((T)0.5*one_over_dy2)/sqrt(POW2(dvpvp) + POW2(delta)))*max_dvv_r;
					u_ul_2 += (T)0.5*((T)0.5*one_over_dy2 - dvpvp*((T)0.5*one_over_dy2)/sqrt(POW2(dvpvp) + POW2(delta)));

					u_ul_3 = 0;
					u_ul_3 = (MA_array[0].array_for_this(i, j) - MA_array[1].array_for_this(i, j))*(u_ul_1 - u_ul_2)/sqrt(POW2(MA_array[0].array_for_this(i, j) - MA_array[1].array_for_this(i, j)) + POW2(delta));
					
					u_ul = 0;
					u_ul = u_ul_1 + u_ul_2 - u_ul_3; 

					// Down-Right Stencil
					T u_dr_1, u_dr_2, u_dr_3;

					u_dr_1 = 0;
										
					u_dr_2 = 0;
					u_dr_2 += (T)0.5*((T)0.5*one_over_dy2 + dvpvp*((T)0.5*one_over_dy2)/sqrt(POW2(dvpvp) + POW2(delta)))*max_dvv_r;
					u_dr_2 += (T)0.5*((T)0.5*one_over_dy2 - dvpvp*((T)0.5*one_over_dy2)/sqrt(POW2(dvpvp) + POW2(delta)));

					u_dr_3 = 0;
					u_dr_3 = (MA_array[0].array_for_this(i, j) - MA_array[1].array_for_this(i, j))*(u_dr_1 - u_dr_2)/sqrt(POW2(MA_array[0].array_for_this(i, j) - MA_array[1].array_for_this(i, j)) + POW2(delta));

					u_dr = 0;
					u_dr = u_dr_1 + u_dr_2 - u_dr_3; 

					if (i - 1 == mid_i && j == mid_j)
					{
						u_l = u_l - 1;
					}
					else if (i + 1 == mid_i && j == mid_j)
					{
						u_r = u_r - 1;
					}
					else if (i == mid_i && j - 1 == mid_j)
					{
						u_d = u_d - 1;
					}
					else if (i == mid_i && j + 1 == mid_j)
					{
						u_u = u_u - 1;
					}
					else if (i == mid_i && j == mid_j)
					{
						u_c = u_c - 1;
					}
					else if (i - 1 == mid_i && j - 1 == mid_j)
					{
						u_dl = u_dl - 1;
					}
					else if (i + 1 == mid_i && j + 1 == mid_j)
					{
						u_ur = u_ur - 1;
					}
					else if (i + 1 == mid_i && j - 1 == mid_j)
					{
						u_dr = u_dr - 1;
					}
					else if (i - 1 == mid_i && j + 1 == mid_j)
					{
						u_ul = u_ul - 1;
					}
					else
					{
						J.AssignValue(bc(i, j), bc(mid_i, mid_j), -1, thread_id);
					}

					if (bc(i - 1, j) > -1)
					{
						J.AssignValue(bc(i, j), bc(i - 1, j), u_l, thread_id);
					}
					if (bc(i + 1, j) > -1)
					{
						J.AssignValue(bc(i, j), bc(i + 1, j), u_r, thread_id);
					}
					if (bc(i, j - 1) > -1)
					{
						J.AssignValue(bc(i, j), bc(i, j - 1), u_d, thread_id);
					}
					if (bc(i, j + 1) > -1)
					{
						J.AssignValue(bc(i, j), bc(i, j + 1), u_u, thread_id);
					}
					
					if (bc(i - 1, j - 1) > -1)
					{
						J.AssignValue(bc(i, j), bc(i - 1, j - 1), u_dl, thread_id);
					}
					if (bc(i + 1, j + 1) > -1)
					{
						J.AssignValue(bc(i, j), bc(i + 1, j + 1), u_ur, thread_id);
					}
					if (bc(i + 1, j - 1) > -1)
					{
						J.AssignValue(bc(i, j), bc(i + 1, j - 1), u_dr, thread_id);
					}
					if (bc(i - 1, j + 1) > -1)
					{
						J.AssignValue(bc(i, j), bc(i - 1, j + 1), u_ul, thread_id);
					}

					J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id); 
				}
				else
				{
					if ((MA_array.values[1])(i, j) >= (MA_array.values[0])(i, j))
					{
						b_vector[bc(i, j)] = MA(i, j);

						T u_l, u_r, u_c, u_d, u_u;

						if (dxx >= delta)
						{
							if (dyy >= delta)
							{
								u_l = one_over_dx2*dyy + partial_p(i, j)*one_over_2dx;
								u_r = one_over_dx2*dyy - partial_p(i, j)*one_over_2dx;
								u_c = -2*one_over_dx2*(dxx + dyy);
								u_d = one_over_dy2*dxx + partial_q(i, j)*one_over_2dy;
								u_u = one_over_dy2*dxx - partial_q(i, j)*one_over_2dy;

								if (i - 1 == mid_i && j == mid_j)
								{
									u_l = u_l - 4;
								}
								else if (i + 1 == mid_i && j == mid_j)
								{
									u_r = u_r - 4;
								}
								else if (i == mid_i && j - 1 == mid_j)
								{
									u_d = u_d - 4;
								}
								else if (i == mid_i && j + 1 == mid_j)
								{
									u_u = u_u - 4;
								}
								else if (i == mid_i && j == mid_j)
								{
									u_c = u_c - 4;
								}
								else
								{
									J.AssignValue(bc(i, j), bc(mid_i, mid_j), -4, thread_id);
								}

								if (bc(i - 1, j) > -1)
								{
									J.AssignValue(bc(i, j), bc(i - 1, j), u_l, thread_id);
								}
								if (bc(i + 1, j) > -1)
								{
									J.AssignValue(bc(i, j), bc(i + 1, j), u_r, thread_id);
								}
								if (bc(i, j - 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i, j - 1), u_d, thread_id);
								}
								if (bc(i, j + 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i, j + 1), u_u, thread_id);
								}

								J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id); 
							}
							else
							{
								u_l = one_over_dx2*delta + partial_p(i, j)*one_over_2dx;
								u_r = one_over_dx2*delta - partial_p(i, j)*one_over_2dx;
								u_c = -2*one_over_dx2*(delta - 1);
								u_d = -one_over_dy2 + partial_q(i, j)*one_over_2dy;
								u_u = -one_over_dy2 - partial_q(i, j)*one_over_2dy;

								if (i - 1 == mid_i && j == mid_j)
								{
									u_l = u_l - 4;
								}
								else if (i + 1 == mid_i && j == mid_j)
								{
									u_r = u_r - 4;
								}
								else if (i == mid_i && j - 1 == mid_j)
								{
									u_d = u_d - 4;
								}
								else if (i == mid_i && j + 1 == mid_j)
								{
									u_u = u_u - 4;
								}
								else if (i == mid_i && j == mid_j)
								{
									u_c = u_c - 4;
								}
								else
								{
									J.AssignValue(bc(i, j), bc(mid_i, mid_j), -4, thread_id);
								}

								if (bc(i - 1, j) > -1)
								{
									J.AssignValue(bc(i, j), bc(i - 1, j), u_l, thread_id);
								}
								if (bc(i + 1, j) > -1)
								{
									J.AssignValue(bc(i, j), bc(i + 1, j), u_r, thread_id);
								}
								if (bc(i, j - 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i, j - 1), u_d, thread_id);
								}
								if (bc(i, j + 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i, j + 1), u_u, thread_id);
								}

								J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
							}
						}
						else
						{
							T u_l, u_r, u_c, u_d, u_u;

							if (dyy >= delta)
							{
								u_l = -one_over_dx2 + partial_p(i, j)*one_over_2dx;
								u_r = -one_over_dx2 - partial_p(i, j)*one_over_2dx;
								u_c = -2*one_over_dx2*(delta - 1);
								u_d = one_over_dy2*delta + partial_q(i, j)*one_over_2dy;
								u_u = one_over_dy2*delta - partial_q(i, j)*one_over_2dy;

								if (i - 1 == mid_i && j == mid_j)
								{
									u_l = u_l - 4;
								}
								else if (i + 1 == mid_i && j == mid_j)
								{
									u_r = u_r - 4;
								}
								else if (i == mid_i && j - 1 == mid_j)
								{
									u_d = u_d - 4;
								}
								else if (i == mid_i && j + 1 == mid_j)
								{
									u_u = u_u - 4;
								}
								else if (i == mid_i && j == mid_j)
								{
									u_c = u_c - 4;
								}
								else
								{
									J.AssignValue(bc(i, j), bc(mid_i, mid_j), -4, thread_id);
								}

								if (bc(i - 1, j) > -1)
								{
									J.AssignValue(bc(i, j), bc(i - 1, j), u_l, thread_id);
								}
								if (bc(i + 1, j) > -1)
								{
									J.AssignValue(bc(i, j), bc(i + 1, j), u_r, thread_id);
								}
								if (bc(i, j - 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i, j - 1), u_d, thread_id);
								}
								if (bc(i, j + 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i, j + 1), u_u, thread_id);
								}

								J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
							}
							else
							{
								u_l = -one_over_dx2 + partial_p(i, j)*one_over_2dx;
								u_r = -one_over_dx2 - partial_p(i, j)*one_over_2dx;
								u_c = 4*one_over_dx2;
								u_d = -one_over_dy2 + partial_q(i, j)*one_over_2dy;
								u_u = -one_over_dy2 - partial_q(i, j)*one_over_2dy;

								if (i - 1 == mid_i && j == mid_j)
								{
									u_l = u_l - 4;
								}
								else if (i + 1 == mid_i && j == mid_j)
								{
									u_r = u_r - 4;
								}
								else if (i == mid_i && j - 1 == mid_j)
								{
									u_d = u_d - 4;
								}
								else if (i == mid_i && j + 1 == mid_j)
								{
									u_u = u_u - 4;
								}
								else if (i == mid_i && j == mid_j)
								{
									u_c = u_c - 1;
								}
								else
								{
									J.AssignValue(bc(i, j), bc(mid_i, mid_j), -4, thread_id);
								}
								
								if (bc(i - 1, j) > -1)
								{
									J.AssignValue(bc(i, j), bc(i - 1, j), u_l, thread_id);
								}
								if (bc(i + 1, j) > -1)
								{
									J.AssignValue(bc(i, j), bc(i + 1, j), u_r, thread_id);
								}
								if (bc(i, j - 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i, j - 1), u_d, thread_id);
								}
								if (bc(i, j + 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i, j + 1), u_u, thread_id);
								}

								J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
							}
						}
					}
					else
					{
						b_vector[bc(i, j)] = MA(i, j);

						T u_dl, u_ur, u_c, u_dr, u_ul;

						if (dvv >= delta)
						{
							if (dvpvp >= delta)
							{
								u_dl = (T)1/2*one_over_dx2*dvpvp + (T)0.5*partial_p(i, j)*one_over_2dx + (T)0.5*partial_q(i, j)*one_over_2dy;
								u_ur = (T)1/2*one_over_dx2*dvpvp - (T)0.5*partial_p(i, j)*one_over_2dx - (T)0.5*partial_q(i, j)*one_over_2dy;
								u_c  = -one_over_dx2*(dvv + dvpvp);
								u_dr = (T)1/2*one_over_dx2*dvv - (T)0.5*partial_p(i, j)*one_over_2dx + (T)0.5*partial_q(i, j)*one_over_2dy;
								u_ul = (T)1/2*one_over_dx2*dvv + (T)0.5*partial_p(i, j)*one_over_2dx - (T)0.5*partial_q(i, j)*one_over_2dy;

								if (i - 1 == mid_i && j - 1 == mid_j)
								{
									u_dl = u_dl - 4;
								}
								else if (i + 1 == mid_i && j + 1 == mid_j)
								{
									u_ur = u_ur - 4;
								}
								else if (i + 1 == mid_i && j - 1 == mid_j)
								{
									u_dr = u_dr - 4;
								}
								else if (i - 1 == mid_i && j + 1 == mid_j)
								{
									u_ul = u_ul - 4;
								}
								else if (i == mid_i && j == mid_j)
								{
									u_c = u_c - 4;
								}
								else
								{
									J.AssignValue(bc(i, j), bc(mid_i, mid_j), -4, thread_id);
								}
								
								if (bc(i - 1, j - 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i - 1, j - 1), u_dl, thread_id);
								}
								if (bc(i + 1, j + 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i + 1, j + 1), u_ur, thread_id);
								}
								if (bc(i + 1, j - 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i + 1, j - 1), u_dr, thread_id);
								}
								if (bc(i - 1, j + 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i - 1, j + 1), u_ul, thread_id);
								}

								J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
							}
							else
							{
								u_dl = (T)1/2*one_over_dx2*delta + (T)0.5*partial_p(i, j)*one_over_2dx + (T)0.5*partial_q(i, j)*one_over_2dy;
								u_ur = (T)1/2*one_over_dx2*delta - (T)0.5*partial_p(i, j)*one_over_2dx - (T)0.5*partial_q(i, j)*one_over_2dy;
								u_c  = -one_over_dx2*(delta - 1);
								u_dr = -(T)1/2*one_over_dx2 - (T)0.5*partial_p(i, j)*one_over_2dx + (T)0.5*partial_q(i, j)*one_over_2dy;
								u_ul = -(T)1/2*one_over_dx2 + (T)0.5*partial_p(i, j)*one_over_2dx - (T)0.5*partial_q(i, j)*one_over_2dy;

								if (i - 1 == mid_i && j - 1 == mid_j)
								{
									u_dl = u_dl - 4;
								}
								else if (i + 1 == mid_i && j + 1 == mid_j)
								{
									u_ur = u_ur - 4;
								}
								else if (i + 1 == mid_i && j - 1 == mid_j)
								{
									u_dr = u_dr - 4;
								}
								else if (i - 1 == mid_i && j + 1 == mid_j)
								{
									u_ul = u_ul - 4;
								}
								else if (i == mid_i && j == mid_j)
								{
									u_c = u_c - 4;
								}
								else
								{
									J.AssignValue(bc(i, j), bc(mid_i, mid_j), -4, thread_id);
								}

								if (bc(i - 1, j - 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i - 1, j - 1), u_dl, thread_id);
								}
								if (bc(i + 1, j + 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i + 1, j + 1), u_ur, thread_id);
								}
								if (bc(i + 1, j - 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i + 1, j - 1), u_dr, thread_id);
								}
								if (bc(i - 1, j + 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i - 1, j + 1), u_ul, thread_id);
								}

								J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
							}
						}
						else
						{
							T u_dl, u_ur, u_c, u_dr, u_ul;

							if (dvpvp >= delta)
							{
								u_dl = -(T)1/2*one_over_dx2 + (T)0.5*partial_p(i, j)*one_over_2dx + (T)0.5*partial_q(i, j)*one_over_2dy;
								u_ur = -(T)1/2*one_over_dx2 - (T)0.5*partial_p(i, j)*one_over_2dx - (T)0.5*partial_q(i, j)*one_over_2dy;
								u_c  = -one_over_dx2*(delta - 1);
								u_dr = (T)1/2*one_over_dx2*delta - (T)0.5*partial_p(i, j)*one_over_2dx + (T)0.5*partial_q(i, j)*one_over_2dy;
								u_ul = (T)1/2*one_over_dx2*delta + (T)0.5*partial_p(i, j)*one_over_2dx - (T)0.5*partial_q(i, j)*one_over_2dy;

								if (i - 1 == mid_i && j - 1 == mid_j)
								{
									u_dl = u_dl - 4;
								}
								else if (i + 1 == mid_i && j + 1 == mid_j)
								{
									u_ur = u_ur - 4;
								}
								else if (i + 1 == mid_i && j - 1 == mid_j)
								{
									u_dr = u_dr - 4;
								}
								else if (i - 1 == mid_i && j + 1 == mid_j)
								{
									u_ul = u_ul - 4;
								}
								else if (i == mid_i && j == mid_j)
								{
									u_c = u_c - 4;
								}
								else
								{
									J.AssignValue(bc(i, j), bc(mid_i, mid_j), -4, thread_id);
								}

								if (bc(i - 1, j - 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i - 1, j - 1), u_dl, thread_id);
								}
								if (bc(i + 1, j + 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i + 1, j + 1), u_ur, thread_id);
								}
								if (bc(i + 1, j - 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i + 1, j - 1), u_dr, thread_id);
								}
								if (bc(i - 1, j + 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i - 1, j + 1), u_ul, thread_id);
								}

								J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
							}
							else
							{
								u_dl = -(T)1/2*one_over_dx2 + (T)0.5*partial_p(i, j)*one_over_2dx + (T)0.5*partial_q(i, j)*one_over_2dy;
								u_ur = -(T)1/2*one_over_dx2 - (T)0.5*partial_p(i, j)*one_over_2dx - (T)0.5*partial_q(i, j)*one_over_2dy;
								u_c  = (T)2*one_over_dx2;
								u_dr = -(T)1/2*one_over_dx2 - (T)0.5*partial_p(i, j)*one_over_2dx + (T)0.5*partial_q(i, j)*one_over_2dy;
								u_ul = -(T)1/2*one_over_dx2 + (T)0.5*partial_p(i, j)*one_over_2dx - (T)0.5*partial_q(i, j)*one_over_2dy;

								if (i - 1 == mid_i && j - 1 == mid_j)
								{
									u_dl = u_dl - 4;
								}
								else if (i + 1 == mid_i && j + 1 == mid_j)
								{
									u_ur = u_ur - 4;
								}
								else if (i + 1 == mid_i && j - 1 == mid_j)
								{
									u_dr = u_dr - 4;
								}
								else if (i - 1 == mid_i && j + 1 == mid_j)
								{
									u_ul = u_ul - 4;
								}
								else if (i == mid_i && j == mid_j)
								{
									u_c = u_c - 4;
								}
								else
								{
									J.AssignValue(bc(i, j), bc(mid_i, mid_j), -4, thread_id);
								}

								if (bc(i - 1, j - 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i - 1, j - 1), u_dl, thread_id);
								}
								if (bc(i + 1, j + 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i + 1, j + 1), u_ur, thread_id);
								}
								if (bc(i + 1, j - 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i + 1, j - 1), u_dr, thread_id);
								}
								if (bc(i - 1, j + 1) > -1)
								{
									J.AssignValue(bc(i, j), bc(i - 1, j + 1), u_ul, thread_id);
								}

								J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
							}
						}
					} 
				} 
			}
		}
		END_GRID_ITERATION_2D;
	}

	int CountNonZeroElements(const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id, const ARRAY<FIELD_STRUCTURE_2D<T>>& stencil)
	{
		int nnz(0);
		int nnz_spe(0);
		int mid_i, mid_j;

		BEGIN_HEAD_THREAD_WORK
		{
			mid_i = (T)0.5*(bc.grid.i_start + bc.grid.i_end);
			mid_j = (T)0.5*(bc.grid.j_start + bc.grid.j_end);
		}
		END_HEAD_THREAD_WORK;

		BEGIN_GRID_ITERATION_2D(bc.partial_grids[thread_id])
		{
			if (bc(i, j) > -1)
			{
				nnz++;
				
				// For the boundary update
				if (i == i_start)
				{
					if (bc(i + 1, j) > -1)
					{
						nnz++;
					}
					if (bc(i - 1, j) > -1)
					{
						nnz++;
					}
					if (bc(i, j + 1) > -1)
					{
						nnz++;
					}
					if (bc(i, j - 1) > -1)
					{
						nnz++;
					} 
				}
				else if (i == i_end)
				{
					if (bc(i + 1, j) > -1)
					{
						nnz++;
					}
					if (bc(i - 1, j) > -1)
					{
						nnz++;
					}
					if (bc(i, j + 1) > -1)
					{
						nnz++;
					}
					if (bc(i, j - 1) > -1)
					{
						nnz++;
					} 
				}
				else if (j == j_start)
				{
					if (bc(i + 1, j) > -1)
					{
						nnz++;
					}
					if (bc(i - 1, j) > -1)
					{
						nnz++;
					}
					if (bc(i, j + 1) > -1)
					{
						nnz++;
					}
					if (bc(i, j - 1) > -1)
					{
						nnz++;
					} 
				}
				else if (j == j_end)
				{
					if (bc(i + 1, j) > -1)
					{
						nnz++;
					}
					if (bc(i - 1, j) > -1)
					{
						nnz++;
					}
					if (bc(i, j + 1) > -1)
					{
						nnz++;
					}
					if (bc(i, j - 1) > -1)
					{
						nnz++;
					} 
				}
				else
				{
					if (is_regularized)
					{
						if (bc(i + 1, j) > -1)
						{
							nnz++;
						}
						if (bc(i - 1, j) > -1)
						{
							nnz++;
						}
						if (bc(i, j + 1) > -1)
						{
							nnz++;
						}
						if (bc(i, j - 1) > -1)
						{
							nnz++;
						}
						if (bc(i - 1, j - 1) > -1)
						{
							nnz++;
						}
						if (bc(i - 1, j + 1) > -1)
						{
							nnz++;
						}
						if (bc(i + 1, j - 1) > -1)
						{
							nnz++;
						}
						if (bc(i + 1, j + 1) > -1)
						{
							nnz++;
						} 
					}
					else
					{
						if (stencil.values[1].array_for_this(i, j) >= stencil.values[0].array_for_this(i, j))
						{
							if (bc(i + 1, j) > -1)
							{
								nnz++;
							}
							if (bc(i - 1, j) > -1)
							{
								nnz++;
							}
							if (bc(i, j + 1) > -1)
							{
								nnz++;
							}
							if (bc(i, j - 1) > -1)
							{
								nnz++;
							} 
						}
						else
						{
							if (bc(i - 1, j - 1) > -1)
							{
								nnz++;
							}
							if (bc(i - 1, j + 1) > -1)
							{
								nnz++;
							}
							if (bc(i + 1, j - 1) > -1)
							{
								nnz++;
							}
							if (bc(i + 1, j + 1) > -1)
							{
								nnz++;
							} 
						} 
					}
				}
			}
		}
		END_GRID_ITERATION_SUM(nnz);
		
		return nnz;
	}

	int CountNonZeroElements(const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id, const ARRAY<FIELD_STRUCTURE_2D<T>>& stencil, const int& N)
	{
		int nnz(0);
		
		AssignStencil(density_x, density_y, density_y, thread_id);
		
		// Speedup Variable
		T one_over_dx  = u.grid.one_over_dx, one_over_dy = u.grid.one_over_dy;
		T one_over_dx2 = u.grid.one_over_dx2, one_over_dy2 = u.grid.one_over_dy2;
		T one_over_2dx = u.grid.one_over_2dx, one_over_2dy = u.grid.one_over_2dy;
		
		T x_min = u.grid.x_min, x_max = u.grid.x_max, dx = u.grid.dx;
		T y_min = u.grid.y_min, y_max = u.grid.y_max, dy = u.grid.dy;

		int i_start = u.grid.i_start, i_end = u.grid.i_end;
		int j_start = u.grid.j_start, j_end = u.grid.j_end;
		
		// For fixed point
		int mid_i(0), mid_j(0);
		
		mid_i = (T)0.5*(i_start + i_end);
		mid_j = (T)0.5*(j_start + j_end);

		T H_star_n(0);

		ARRAY<VT> target_normal;
		target_normal.Initialize(n_y);

		for (int i = 0; i < n_y; i++)
		{
			target_normal[i].x = cos(2 * PI*i / n_y);
			target_normal[i].y = sin(2 * PI*i / n_y);
		}

		// Define H*(n)
		H_n.Initialize(n_y);

		if (test_number == 1)
		{
			for (int k = 0; k < target_normal.length; k++)
			{
				VT& tn(target_normal[k]);

				T x_coor_0 = x_min + i_start*dx, y_coor_0 = y_min + j_start*dy;

				VT y_start(x_coor_0, y_coor_0);

				H_n[k] = DotProduct(y_start, tn);

				// Left Boundary
				for (int j = j_start; j <= j_end; j++)
				{
					T x_coor = x_min + i_start*dx, y_coor = y_min + j*dy;

					VT y_0(x_coor, y_coor);

					H_n[k] = MAX(H_n[k], DotProduct(y_0, tn));
				}

				// Right Boundary
				for (int j = j_start; j <= j_end; j++)
				{
					T x_coor = x_min + i_end*dx, y_coor = y_min + j*dy;

					VT y_0(x_coor, y_coor);

					H_n[k] = MAX(H_n[k], DotProduct(y_0, tn));
				}

				// Lower Boundary
				for (int i = i_start; i <= i_end; i++)
				{
					T x_coor = x_min + i*dx, y_coor = y_min + j_start*dy;

					VT y_0(x_coor, y_coor);

					H_n[k] = MAX(H_n[k], DotProduct(y_0, tn));
				}

				// Upper Boundary
				for (int i = i_start; i <= i_end; i++)
				{
					T x_coor = x_min + i*dx, y_coor = y_min + j_end*dy;

					VT y_0(x_coor, y_coor);

					H_n[k] = MAX(H_n[k], DotProduct(y_0, tn));
				} 
			} 
		}

		if (test_number == 2)
		{
			/*int num_of_target_boundary(0);

			BEGIN_GRID_ITERATION_2D(u.partial_grids[thread_id])
			{
				T x_coor = u.grid.x_min + i*u.grid.dx, y_coor = u.grid.y_min + j*u.grid.dy;
				T x_coor_l = u.grid.x_min + (i - 1)*u.grid.dx, y_coor_d = u.grid.y_min + (j - 1)*u.grid.dy;
				T x_coor_r = u.grid.x_min + (i + 1)*u.grid.dx, y_coor_u = u.grid.y_min + (j + 1)*u.grid.dy;
				
				T domain_region = POW2(x_coor)/POW2(0.8) + POW2(y_coor)/POW2(0.4) - 1;
				T domain_region_r = POW2(x_coor_r)/POW2(0.8) + POW2(y_coor)/POW2(0.4) - 1;
				T domain_region_l = POW2(x_coor_l)/POW2(0.8) + POW2(y_coor)/POW2(0.4) - 1;
				T domain_region_u = POW2(x_coor)/POW2(0.8) + POW2(y_coor_u)/POW2(0.4) - 1;
				T domain_region_d = POW2(x_coor)/POW2(0.8) + POW2(y_coor_d)/POW2(0.4) - 1;

				if (domain_region < 0)
				{
					if (domain_region_r > 0)
					{
						num_of_target_boundary++;
					}
					else if (domain_region_l > 0)
					{
						num_of_target_boundary++;
					}
					else if (domain_region_d > 0)
					{
						num_of_target_boundary++;
					}
					else if (domain_region_u > 0)
					{
						num_of_target_boundary++;
					}
				}
				if (domain_region == 0)
				{
					num_of_target_boundary++;
				}
			}
			END_GRID_ITERATION_2D;
			
			ARRAY<VT> target_boundary;
			target_boundary.Initialize(num_of_target_boundary);
			
			int index_for_target_boundary(0);

			BEGIN_GRID_ITERATION_2D(u.partial_grids[thread_id])
			{
				T x_coor = u.grid.x_min + i*u.grid.dx, y_coor = u.grid.y_min + j*u.grid.dy;
				T x_coor_l = u.grid.x_min + (i - 1)*u.grid.dx, y_coor_d = u.grid.y_min + (j - 1)*u.grid.dy;
				T x_coor_r = u.grid.x_min + (i + 1)*u.grid.dx, y_coor_u = u.grid.y_min + (j + 1)*u.grid.dy;

				T domain_region = POW2(x_coor)/POW2(0.8) + POW2(y_coor)/POW2(0.4) - 1;
				T domain_region_r = POW2(x_coor_r)/POW2(0.8) + POW2(y_coor)/POW2(0.4) - 1;
				T domain_region_l = POW2(x_coor_l)/POW2(0.8) + POW2(y_coor)/POW2(0.4) - 1;
				T domain_region_u = POW2(x_coor)/POW2(0.8) + POW2(y_coor_u)/POW2(0.4) - 1;
				T domain_region_d = POW2(x_coor)/POW2(0.8) + POW2(y_coor_d)/POW2(0.4) - 1;

				MATRIX_2X2 M_x(0.8, 0, 0, 0.4), M_y(0.6, 0.2, 0.2, 0.8), J(0, 1, -1, 0);

				MATRIX_2X2 nu, deno;

				nu = M_x.Inversed()*M_y.Inversed()*J;
				deno = M_x.Inversed()*M_y.Inversed();
			
				T atheta = nu.Trace()/deno.Trace();

				MATRIX_2X2 R_theta(cos(atan(atheta)), sin(atan(atheta)), -sin(atan(atheta)), cos(atan(atheta)));

				MATRIX_2X2 Ellipse_Mapping;

				Ellipse_Mapping = M_y*R_theta*M_x.Inversed();

				T x_vec_c = Ellipse_Mapping.x[0]*x_coor + Ellipse_Mapping.x[2]*y_coor, y_vec_c = Ellipse_Mapping.x[1]*x_coor + Ellipse_Mapping.x[3]*y_coor;

				if (domain_region < 0)
				{
					if (domain_region_r > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_vec_c, y_vec_c);
						index_for_target_boundary++;
					}
					else if (domain_region_l > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_vec_c, y_vec_c);
						index_for_target_boundary++;
					}
					else if (domain_region_d > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_vec_c, y_vec_c);
						index_for_target_boundary++;
					}
					else if (domain_region_u > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_vec_c, y_vec_c);
						index_for_target_boundary++;
					}
				}
				if (domain_region == 0)
				{
					target_boundary[index_for_target_boundary] = VT(x_vec_c, y_vec_c);
					index_for_target_boundary++;
				}
			}
			END_GRID_ITERATION_2D;

			for (int k = 0; k < n_y; k++)
			{
				VT& tn(target_normal[k]);

				H_n[k] = DotProduct(target_boundary[0], tn);
				
				for (int i = 0; i < num_of_target_boundary; i++)
				{
					H_n[k] = MAX(H_n[k], DotProduct(target_boundary[i], tn));
				}
			}*/
			
			// Boundary Vector define
			MATRIX_2X2 M_x(0.8, 0, 0, 0.4), M_y(0.6, 0.2, 0.2, 0.8), J(0, 1, -1, 0);

			MATRIX_2X2 nu, deno;

			nu = M_x.Inversed()*M_y.Inversed()*J;
			deno = M_x.Inversed()*M_y.Inversed();
			
			T atheta = nu.Trace()/deno.Trace();

			MATRIX_2X2 R_theta(cos(atan(atheta)), sin(atan(atheta)), -sin(atan(atheta)), cos(atan(atheta)));

			MATRIX_2X2 Ellipse_Mapping;

			Ellipse_Mapping = M_y*R_theta*M_x.Inversed();

			ARRAY<VT> target_boundary;
			target_boundary.Initialize(N_point);
				
			for (int i = 0; i < N_point; i++)
			{
				T angle = 2*i*PI/N_point;

				target_boundary[i].x = 0.8*Ellipse_Mapping.x[0]*cos(angle) + 0.4*Ellipse_Mapping.x[2]*sin(angle);
				target_boundary[i].y = 0.8*Ellipse_Mapping.x[1]*cos(angle) + 0.4*Ellipse_Mapping.x[3]*sin(angle);

				/*target_boundary[i].x = 0.6*cos(angle) - 0.2*sin(angle);
				target_boundary[i].y = 0.2*cos(angle) - 0.8*sin(angle);*/
			}
				
			for (int k = 0; k < n_y; k++)
			{
				VT& tn(target_normal[k]);

				H_n[k] = DotProduct(target_boundary[0], tn);
				
				for (int i = 0; i < N_point; i++)
				{
					H_n[k] = MAX(H_n[k], DotProduct(target_boundary[i], tn));
				}
			}
		}

		if (test_number == 3)
		{
			int num_of_target_boundary(0);

			BEGIN_GRID_ITERATION_2D(u.partial_grids[thread_id])
			{
				T x_coor = u.grid.x_min + i*u.grid.dx, y_coor = u.grid.y_min + j*u.grid.dy;
				T x_coor_l = u.grid.x_min + (i - 1)*u.grid.dx, y_coor_d = u.grid.y_min + (j - 1)*u.grid.dy;
				T x_coor_r = u.grid.x_min + (i + 1)*u.grid.dx, y_coor_u = u.grid.y_min + (j + 1)*u.grid.dy;

				T domain_region_1 = POW2(x_coor + 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_1_r = POW2(x_coor_r + 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_1_l = POW2(x_coor_l + 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_1_u = POW2(x_coor + 0.1) + POW2(y_coor_u) - POW2(0.85);
				T domain_region_1_d = POW2(x_coor + 0.1) + POW2(y_coor_d) - POW2(0.85);

				T domain_region_2 = POW2(x_coor - 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_2_r = POW2(x_coor_r - 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_2_l = POW2(x_coor_l - 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_2_u = POW2(x_coor - 0.1) + POW2(y_coor_u) - POW2(0.85);
				T domain_region_2_d = POW2(x_coor - 0.1) + POW2(y_coor_d) - POW2(0.85);

				if (domain_region_1 < 0 && x_coor < -0.1)
				{
					if (domain_region_1_r > 0)
					{
						num_of_target_boundary++;
					}
					else if (domain_region_1_l > 0)
					{
						num_of_target_boundary++;
					}
					else if (domain_region_1_d > 0)
					{
						num_of_target_boundary++;
					}
					else if (domain_region_1_u > 0)
					{
						num_of_target_boundary++;
					}
				}
				if (domain_region_2 < 0 && x_coor > 0.1)
				{
					if (domain_region_2_r > 0)
					{
						num_of_target_boundary++;
					}
					else if (domain_region_2_l > 0)
					{
						num_of_target_boundary++;
					}
					else if (domain_region_2_d > 0)
					{
						num_of_target_boundary++;
					}
					else if (domain_region_2_u > 0)
					{
						num_of_target_boundary++;
					}
				}
			}
			END_GRID_ITERATION_2D;
			
			ARRAY<VT> target_boundary;
			target_boundary.Initialize(num_of_target_boundary);
			
			int index_for_target_boundary(0);

			BEGIN_GRID_ITERATION_2D(u.partial_grids[thread_id])
			{
				T x_coor = u.grid.x_min + i*u.grid.dx, y_coor = u.grid.y_min + j*u.grid.dy;
				T x_coor_l = u.grid.x_min + (i - 1)*u.grid.dx, y_coor_d = u.grid.y_min + (j - 1)*u.grid.dy;
				T x_coor_r = u.grid.x_min + (i + 1)*u.grid.dx, y_coor_u = u.grid.y_min + (j + 1)*u.grid.dy;

				T domain_region_1 = POW2(x_coor + 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_1_r = POW2(x_coor_r + 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_1_l = POW2(x_coor_l + 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_1_u = POW2(x_coor + 0.1) + POW2(y_coor_u) - POW2(0.85);
				T domain_region_1_d = POW2(x_coor + 0.1) + POW2(y_coor_d) - POW2(0.85);

				T domain_region_2 = POW2(x_coor - 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_2_r = POW2(x_coor_r - 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_2_l = POW2(x_coor_l - 0.1) + POW2(y_coor) - POW2(0.85);
				T domain_region_2_u = POW2(x_coor - 0.1) + POW2(y_coor_u) - POW2(0.85);
				T domain_region_2_d = POW2(x_coor - 0.1) + POW2(y_coor_d) - POW2(0.85);

				if (domain_region_1 < 0 && x_coor < -0.1)
				{
					if (domain_region_1_r > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_coor + 0.1, y_coor);
						index_for_target_boundary++;
					}
					else if (domain_region_1_l > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_coor + 0.1, y_coor);
						index_for_target_boundary++;
					}
					else if (domain_region_1_d > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_coor + 0.1, y_coor);
						index_for_target_boundary++;
					}
					else if (domain_region_1_u > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_coor + 0.1, y_coor);
						index_for_target_boundary++;
					}
				}
				if (domain_region_2 < 0 && x_coor > 0.1)
				{
					if (domain_region_2_r > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_coor - 0.1, y_coor);
						index_for_target_boundary++;
					}
					else if (domain_region_2_l > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_coor - 0.1, y_coor);
						index_for_target_boundary++;
					}
					else if (domain_region_2_d > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_coor - 0.1, y_coor);
						index_for_target_boundary++;
					}
					else if (domain_region_2_u > 0)
					{
						target_boundary[index_for_target_boundary] = VT(x_coor - 0.1, y_coor);
						index_for_target_boundary++;
					}
				}
			}
			END_GRID_ITERATION_2D;

			for (int k = 0; k < n_y; k++)
			{
				VT& tn(target_normal[k]);

				H_n[k] = DotProduct(target_boundary[0], tn);
				
				for (int i = 0; i < num_of_target_boundary; i++)
				{
					H_n[k] = MAX(H_n[k], DotProduct(target_boundary[i], tn));
				}
			}

			// Boundary Vector define
			/*ARRAY<VT> target_boundary;
			target_boundary.Initialize(N_point);
				
			for (int i = 0; i < N_point; i++)
			{
				T angle = 2*i*PI/N_point;

				target_boundary[i].x = 0.85*cos(angle);
				target_boundary[i].y = 0.85*sin(angle);
			}
				
			for (int k = 0; k < n_y; k++)
			{
				VT& tn(target_normal[k]);

				H_n[k] = DotProduct(target_boundary[0], tn);
				
				for (int i = 0; i < N_point; i++)
				{
					H_n[k] = MAX(H_n[k], DotProduct(target_boundary[i], tn));
				}
			}*/
		}

		if (test_number == 8)
		{
			// Boundary Vector define
			T a = cos(PI/8), radius = sqrt((T)1/POW2(a) - 1);
			
			ARRAY<VT> target_boundary;
			target_boundary.Initialize(N_point);
				
			for (int i = 0; i < N_point; i++)
			{
				target_boundary[i].x = radius*cos(2*i*PI/N_point);
				target_boundary[i].y = radius*sin(2*i*PI/N_point) + (T)1/a;
			}
				
			for (int k = 0; k < n_y; k++)
			{
				VT& tn(target_normal[k]);

				H_n[k] = DotProduct(target_boundary[0], tn);
				
				for (int i = 0; i < N_point; i++)
				{
					H_n[k] = MAX(H_n[k], DotProduct(target_boundary[i], tn));
				}
			}
		}
				
		BEGIN_GRID_ITERATION_2D(u.partial_grids[thread_id])
		{
			T dxx, dyy, dvv, dvpvp;

			dxx = one_over_dx2*(u(i + 1, j) - 2*u(i, j) + u(i - 1, j));
			dyy = one_over_dy2*(u(i, j + 1) - 2*u(i, j) + u(i, j - 1));
			dvv   = (T)1/2*one_over_dx2*(u(i + 1, j + 1) + u(i - 1, j - 1) - 2*u(i, j));
			dvpvp = (T)1/2*one_over_dx2*(u(i + 1, j - 1) + u(i - 1, j + 1) - 2*u(i, j));

			if (bc(i, j) < 0)
			{
				continue;
			}
						
			// Need to add the derivative of source term
			if (bc(i - 1, j) < 0)
			{
				if (bc(i, j - 1) < 0)
				{
					if (bc(i + 1, j) > -1)
					{
						nnz++;
					}
					if (bc(i, j + 1) > -1)
					{
						nnz++;
					}	
					nnz++;
				}
				else if (bc(i, j + 1) < 0)
				{
					if (bc(i + 1, j) > -1)
					{
						nnz++;
					}
					if (bc(i, j - 1) > -1)
					{
						nnz++;
					}	

					nnz++;
				}
				else
				{
					// Define b_vector
					T dxp = one_over_dx*(u(i + 1, j) - u(i, j)), dym = one_over_dy*(u(i, j) - u(i, j - 1)), dyp = one_over_dy*(u(i, j + 1) - u(i, j));

					T b_vector_value(0);
					int sub_k(0);
					int max_left(0);

					for (int k = 0; k < target_normal.length; k++)
					{
						VT n_x(-1, 0);
						VT& tn(target_normal[k]);

						T dot = DotProduct(n_x, tn);

						if (dot <= 0)
						{
							continue;
						}
						else
						{
							if (sub_k == 0)
							{
								b_vector_value = tn.x*dxp + MAX(tn.y, 0)*dym + MIN(tn.y, 0)*dyp - H_n[k];
							}

							T compared_one = tn.x*dxp + MAX(tn.y, 0)*dym + MIN(tn.y, 0)*dyp - H_n[k];
						
							if (b_vector_value <= compared_one)
							{
								b_vector_value = compared_one;
								max_left = k;
							}
						
							sub_k++;
						}
					}

					if (target_normal[max_left].y > 0)
					{
						if (bc(i + 1, j) > -1)
						{
							nnz++;
						}
						if (bc(i, j - 1) > -1)
						{
							nnz++;
						}
						
						nnz++;
					}
					else if (target_normal[max_left].y < 0)
					{
						if (bc(i + 1, j) > -1)
						{
							nnz++;
						}
						if (bc(i, j + 1) > -1)
						{
							nnz++;
						}

						nnz++;
					}
					else
					{
						if (bc(i + 1, j) > -1)
						{
							nnz++;
						}
						
						nnz++;
					}
				}
			}
			else if (bc(i + 1, j) < 0)
			{
				if (bc(i, j - 1) < 0)
				{
					if (bc(i - 1, j) > -1)
					{
						nnz++;
					}
					if (bc(i, j + 1) > -1)
					{
						nnz++;
					}	
				
					nnz++;
				}
				else if (bc(i, j + 1) < 0)
				{
					if (bc(i - 1, j) > -1)
					{
						nnz++;
					}
					if (bc(i, j - 1) > -1)
					{
						nnz++;
					}	
				
					nnz++;
				}
				else
				{
					// Define b_vector
					T dxm = one_over_dx*(u(i, j) - u(i - 1, j)), dym = one_over_dy*(u(i, j) - u(i, j - 1)), dyp = one_over_dy*(u(i, j + 1) - u(i, j));

					T b_vector_value(0);
					int max_right(0);
					int sub_k(0);

					for (int k = 0; k < target_normal.length; k++)
					{
						VT n_x(1, 0);
						VT& tn(target_normal[k]);

						T dot = DotProduct(n_x, tn);

						if (dot <= 0)
						{
							continue;
						}

						if (sub_k == 0)
						{
							b_vector_value = tn.x*dxm + MAX(tn.y, 0)*dym + MIN(tn.y, 0)*dyp - H_n[k];
						}

						T compared_one = tn.x*dxm + MAX(tn.y, 0)*dym + MIN(tn.y, 0)*dyp - H_n[k];

						if (b_vector_value <= compared_one)
						{
							b_vector_value = compared_one;
							max_right = k;
						}
						sub_k++;
					}

					if (target_normal[max_right].y > 0)
					{
						if (bc(i - 1, j) > -1)
						{
							nnz++;
						}
						if (bc(i, j - 1) > -1)
						{
							nnz++;
						}
						
						nnz++;
					}
					else if (target_normal[max_right].y < 0)
					{
						if (bc(i - 1, j) > -1)
						{
							nnz++;
						}
						if (bc(i, j + 1) > -1)
						{
							nnz++;
						}

						nnz++;
					}
					else
					{
						if (bc(i - 1, j) > -1)
						{
							nnz++;
						}
						
						nnz++;
					}
				}
			}
			else if (bc(i, j - 1) < 0)
			{
				// Define b_vector
				T dxm = one_over_dx*(u(i, j) - u(i - 1, j)), dxp = one_over_dx*(u(i + 1, j) - u(i, j)), dyp = one_over_dy*(u(i, j + 1) - u(i, j));

				T b_vector_value(0);
				int sub_k(0);

				for (int k = 0; k < target_normal.length; k++)
				{
					VT n_y(0, -1);
					VT& tn(target_normal[k]);

					T dot = DotProduct(n_y, tn);

					if (dot <= 0)
					{
						continue;
					}

					if (sub_k == 0)
					{
						b_vector_value = MAX(tn.x, 0)*dxm + MIN(tn.x, 0)*dxp + tn.y*dyp - H_n[k];
					}

					T compared_one = MAX(tn.x, 0)*dxm + MIN(tn.x, 0)*dxp + tn.y*dyp - H_n[k];

					if (b_vector_value <= compared_one)
					{
						b_vector_value = compared_one;
						max_lower = k;
					}
					sub_k++;
				}

				if (target_normal[max_lower].x > 0)
				{
					if (bc(i - 1, j) > -1)
					{
						nnz++;
					}
					if (bc(i, j + 1) > -1)
					{
						nnz++;
					}
						
					nnz++;
				}
				else if (target_normal[max_lower].x < 0)
				{
					if (bc(i + 1, j) > -1)
					{
						nnz++;
					}
					if (bc(i, j + 1) > -1)
					{
						nnz++;
					}

					nnz++;
				}
				else
				{
					if (bc(i, j + 1) > -1)
					{
						nnz++;
					}
						
					nnz++;
				}
			}
			else if (bc(i, j + 1) < 0)
			{
				// Define b_vector
				T dxm = one_over_dx*(u(i, j) - u(i - 1, j)), dxp = one_over_dx*(u(i + 1, j) - u(i, j)), dym = one_over_dy*(u(i, j) - u(i, j - 1));

				T b_vector_value = 0;
				int sub_k(0);

				for (int k = 0; k < target_normal.length; k++)
				{
					VT n_y(0, 1);
					VT& tn(target_normal[k]);

					T dot = DotProduct(n_y, tn);

					if (dot <= 0)
					{
						continue;
					}

					if (sub_k == 0)
					{
						b_vector_value = MAX(tn.x, 0)*dxm + MIN(tn.x, 0)*dxp + tn.y*dym - H_n[k];
					}

					T compared_one = MAX(tn.x, 0)*dxm + MIN(tn.x, 0)*dxp + tn.y*dym - H_n[k];

					if (b_vector_value <= compared_one)
					{
						b_vector_value = compared_one;
						max_upper = k;
					}
					sub_k++;
				}

				if (target_normal[max_upper].x > 0)
				{
					if (bc(i - 1, j) > -1)
					{
						nnz++;
					}
					if (bc(i, j - 1) > -1)
					{
						nnz++;
					}
						
					nnz++;
				}
				else if (target_normal[max_upper].x < 0)
				{
					if (bc(i + 1, j) > -1)
					{
						nnz++;
					}
					if (bc(i, j - 1) > -1)
					{
						nnz++;
					}

					nnz++;
				}
				else
				{
					if (bc(i, j - 1) > -1)
					{
						nnz++;
					}
						
					nnz++;
				}
			}
			else
			{
				if (is_regularized)
				{
					if (i - 1 == mid_i && j == mid_j)
					{
						
					}
					else if (i + 1 == mid_i && j == mid_j)
					{
						
					}
					else if (i == mid_i && j - 1 == mid_j)
					{
						
					}
					else if (i == mid_i && j + 1 == mid_j)
					{
						
					}
					else if (i == mid_i && j == mid_j)
					{
						
					}
					else if (i - 1 == mid_i && j - 1 == mid_j)
					{
						
					}
					else if (i + 1 == mid_i && j + 1 == mid_j)
					{
						
					}
					else if (i + 1 == mid_i && j - 1 == mid_j)
					{
						
					}
					else if (i - 1 == mid_i && j + 1 == mid_j)
					{
						
					}
					else if (i == mid_i && j == mid_j)
					{
						
					}
					else
					{
						nnz++;
					}

					if (bc(i - 1, j) > -1)
					{
						nnz++;
					}
					if (bc(i + 1, j) > -1)
					{
						nnz++;
					}
					if (bc(i, j - 1) > -1)
					{
						nnz++;
					}
					if (bc(i, j + 1) > -1)
					{
						nnz++;
					}
					
					if (bc(i - 1, j - 1) > -1)
					{
						nnz++;
					}
					if (bc(i + 1, j + 1) > -1)
					{
						nnz++;
					}
					if (bc(i + 1, j - 1) > -1)
					{
						nnz++;
					}
					if (bc(i - 1, j + 1) > -1)
					{
						nnz++;
					}

					nnz++;
				}
				else
				{
					if ((MA_array.values[1])(i, j) >= (MA_array.values[0])(i, j))
					{
						if (dxx >= delta)
						{
							if (dyy >= delta)
							{
								if (i - 1 == mid_i && j == mid_j)
								{
									
								}
								else if (i + 1 == mid_i && j == mid_j)
								{
									
								}
								else if (i == mid_i && j - 1 == mid_j)
								{
									
								}
								else if (i == mid_i && j + 1 == mid_j)
								{
									
								}
								else if (i == mid_i && j == mid_j)
								{
									
								}
								else
								{
									nnz++;
								}

								if (bc(i - 1, j) > -1)
								{
									nnz++;
								}
								if (bc(i + 1, j) > -1)
								{
									nnz++;
								}
								if (bc(i, j - 1) > -1)
								{
									nnz++;
								}
								if (bc(i, j + 1) > -1)
								{
									nnz++;
								}

								nnz++;
							}
							else
							{
								if (i - 1 == mid_i && j == mid_j)
								{
									
								}
								else if (i + 1 == mid_i && j == mid_j)
								{
									
								}
								else if (i == mid_i && j - 1 == mid_j)
								{
									
								}
								else if (i == mid_i && j + 1 == mid_j)
								{
									
								}
								else if (i == mid_i && j == mid_j)
								{
									
								}
								else
								{
									nnz++;
								}

								if (bc(i - 1, j) > -1)
								{
									nnz++;
								}
								if (bc(i + 1, j) > -1)
								{
									nnz++;
								}
								if (bc(i, j - 1) > -1)
								{
									nnz++;
								}
								if (bc(i, j + 1) > -1)
								{
									nnz++;
								}

								nnz++;
							}
						}
						else
						{
							if (dyy >= delta)
							{
								if (i - 1 == mid_i && j == mid_j)
								{
									
								}
								else if (i + 1 == mid_i && j == mid_j)
								{
									
								}
								else if (i == mid_i && j - 1 == mid_j)
								{
									
								}
								else if (i == mid_i && j + 1 == mid_j)
								{
								
								}
								else if (i == mid_i && j == mid_j)
								{
									
								}
								else
								{
									nnz++;
								}

								if (bc(i - 1, j) > -1)
								{
									nnz++;
								}
								if (bc(i + 1, j) > -1)
								{
									nnz++;
								}
								if (bc(i, j - 1) > -1)
								{
									nnz++;
								}
								if (bc(i, j + 1) > -1)
								{
									nnz++;
								}

								nnz++;
							}
							else
							{
								if (i - 1 == mid_i && j == mid_j)
								{
								}
								else if (i + 1 == mid_i && j == mid_j)
								{
								}
								else if (i == mid_i && j - 1 == mid_j)
								{
								}
								else if (i == mid_i && j + 1 == mid_j)
								{
								}
								else if (i == mid_i && j == mid_j)
								{
								}
								else
								{
									nnz++;
								}
								
								if (bc(i - 1, j) > -1)
								{
									nnz++;
								}
								if (bc(i + 1, j) > -1)
								{
									nnz++;
								}
								if (bc(i, j - 1) > -1)
								{
									nnz++;
								}
								if (bc(i, j + 1) > -1)
								{
									nnz++;
								}

								nnz++;
							}
						}
					}
					else
					{
						if (dvv >= delta)
						{
							if (dvpvp >= delta)
							{
								if (i - 1 == mid_i && j - 1 == mid_j)
								{
									
								}
								else if (i + 1 == mid_i && j + 1 == mid_j)
								{
									
								}
								else if (i + 1 == mid_i && j - 1 == mid_j)
								{
									
								}
								else if (i - 1 == mid_i && j + 1 == mid_j)
								{
									
								}
								else if (i == mid_i && j == mid_j)
								{
									
								}
								else
								{
									nnz++;
								}
								
								if (bc(i - 1, j - 1) > -1)
								{
									nnz++;
								}
								if (bc(i + 1, j + 1) > -1)
								{
									nnz++;
								}
								if (bc(i + 1, j - 1) > -1)
								{
									nnz++;
								}
								if (bc(i - 1, j + 1) > -1)
								{
									nnz++;
								}

								nnz++;
							}
							else
							{
								if (i - 1 == mid_i && j - 1 == mid_j)
								{
									
								}
								else if (i + 1 == mid_i && j + 1 == mid_j)
								{
									
								}
								else if (i + 1 == mid_i && j - 1 == mid_j)
								{
									
								}
								else if (i - 1 == mid_i && j + 1 == mid_j)
								{
									
								}
								else if (i == mid_i && j == mid_j)
								{
									
								}
								else
								{
									nnz++;
								}
								
								if (bc(i - 1, j - 1) > -1)
								{
									nnz++;
								}
								if (bc(i + 1, j + 1) > -1)
								{
									nnz++;
								}
								if (bc(i + 1, j - 1) > -1)
								{
									nnz++;
								}
								if (bc(i - 1, j + 1) > -1)
								{
									nnz++;
								}

								nnz++;
							}
						}
						else
						{
							if (dvpvp >= delta)
							{
								if (i - 1 == mid_i && j - 1 == mid_j)
								{
									
								}
								else if (i + 1 == mid_i && j + 1 == mid_j)
								{
									
								}
								else if (i + 1 == mid_i && j - 1 == mid_j)
								{
									
								}
								else if (i - 1 == mid_i && j + 1 == mid_j)
								{
									
								}
								else if (i == mid_i && j == mid_j)
								{
									
								}
								else
								{
									nnz++;
								}
								
								if (bc(i - 1, j - 1) > -1)
								{
									nnz++;
								}
								if (bc(i + 1, j + 1) > -1)
								{
									nnz++;
								}
								if (bc(i + 1, j - 1) > -1)
								{
									nnz++;
								}
								if (bc(i - 1, j + 1) > -1)
								{
									nnz++;
								}

								nnz++;
							}
							else
							{
								if (i - 1 == mid_i && j - 1 == mid_j)
								{
									
								}
								else if (i + 1 == mid_i && j + 1 == mid_j)
								{
									
								}
								else if (i + 1 == mid_i && j - 1 == mid_j)
								{
									
								}
								else if (i - 1 == mid_i && j + 1 == mid_j)
								{
									
								}
								else if (i == mid_i && j == mid_j)
								{
									
								}
								else
								{
									nnz++;
								}
								
								if (bc(i - 1, j - 1) > -1)
								{
									nnz++;
								}
								if (bc(i + 1, j + 1) > -1)
								{
									nnz++;
								}
								if (bc(i + 1, j - 1) > -1)
								{
									nnz++;
								}
								if (bc(i - 1, j + 1) > -1)
								{
									nnz++;
								}

								nnz++;
							}
						}
					} 
				} 
			}
		}
		END_GRID_ITERATION_SUM(nnz);

		return nnz;
	}

	//int CountNonZeroElements(const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id, const ARRAY<FIELD_STRUCTURE_2D<T>>& stencil, const int& N)
	//{
	//	int nnz(0);
	//	int nnz_spe(0);
	//	int mid_i, mid_j;

	//	BEGIN_HEAD_THREAD_WORK
	//	{
	//		mid_i = (T)0.5*(bc.grid.i_start + bc.grid.i_end);
	//		mid_j = (T)0.5*(bc.grid.j_start + bc.grid.j_end);
	//	}
	//	END_HEAD_THREAD_WORK;

	//	BEGIN_GRID_ITERATION_2D(bc.partial_grids[thread_id])
	//	{
	//		if (bc(i, j) > -1)
	//		{
	//			nnz++;
	//			
	//			// For the boundary update
	//			if (i == i_start)
	//			{
	//				if (bc(i + 1, j) > -1)
	//				{
	//					nnz++;
	//				}
	//				if (bc(i - 1, j) > -1)
	//				{
	//					nnz++;
	//				}
	//				if (bc(i, j + 1) > -1)
	//				{
	//					nnz++;
	//				}
	//				if (bc(i, j - 1) > -1)
	//				{
	//					nnz++;
	//				} 
	//			}
	//			else if (i == i_end)
	//			{
	//				if (bc(i + 1, j) > -1)
	//				{
	//					nnz++;
	//				}
	//				if (bc(i - 1, j) > -1)
	//				{
	//					nnz++;
	//				}
	//				if (bc(i, j + 1) > -1)
	//				{
	//					nnz++;
	//				}
	//				if (bc(i, j - 1) > -1)
	//				{
	//					nnz++;
	//				} 
	//			}
	//			else if (j == j_start)
	//			{
	//				if (bc(i + 1, j) > -1)
	//				{
	//					nnz++;
	//				}
	//				if (bc(i - 1, j) > -1)
	//				{
	//					nnz++;
	//				}
	//				if (bc(i, j + 1) > -1)
	//				{
	//					nnz++;
	//				}
	//				if (bc(i, j - 1) > -1)
	//				{
	//					nnz++;
	//				} 
	//			}
	//			else if (j == j_end)
	//			{
	//				if (bc(i + 1, j) > -1)
	//				{
	//					nnz++;
	//				}
	//				if (bc(i - 1, j) > -1)
	//				{
	//					nnz++;
	//				}
	//				if (bc(i, j + 1) > -1)
	//				{
	//					nnz++;
	//				}
	//				if (bc(i, j - 1) > -1)
	//				{
	//					nnz++;
	//				} 
	//			}
	//			else
	//			{
	//				if (is_regularized)
	//				{
	//					if (bc(i + 1, j) > -1)
	//					{
	//						nnz++;
	//					}
	//					if (bc(i - 1, j) > -1)
	//					{
	//						nnz++;
	//					}
	//					if (bc(i, j + 1) > -1)
	//					{
	//						nnz++;
	//					}
	//					if (bc(i, j - 1) > -1)
	//					{
	//						nnz++;
	//					}
	//					if (bc(i - 1, j - 1) > -1)
	//					{
	//						nnz++;
	//					}
	//					if (bc(i - 1, j + 1) > -1)
	//					{
	//						nnz++;
	//					}
	//					if (bc(i + 1, j - 1) > -1)
	//					{
	//						nnz++;
	//					}
	//					if (bc(i + 1, j + 1) > -1)
	//					{
	//						nnz++;
	//					} 
	//				}
	//				else
	//				{
	//					if (stencil.values[1].array_for_this(i, j) >= stencil.values[0].array_for_this(i, j))
	//					{
	//						if (bc(i + 1, j) > -1)
	//						{
	//							nnz++;
	//						}
	//						if (bc(i - 1, j) > -1)
	//						{
	//							nnz++;
	//						}
	//						if (bc(i, j + 1) > -1)
	//						{
	//							nnz++;
	//						}
	//						if (bc(i, j - 1) > -1)
	//						{
	//							nnz++;
	//						} 
	//					}
	//					else
	//					{
	//						if (bc(i - 1, j - 1) > -1)
	//						{
	//							nnz++;
	//						}
	//						if (bc(i - 1, j + 1) > -1)
	//						{
	//							nnz++;
	//						}
	//						if (bc(i + 1, j - 1) > -1)
	//						{
	//							nnz++;
	//						}
	//						if (bc(i + 1, j + 1) > -1)
	//						{
	//							nnz++;
	//						} 
	//					} 
	//				}
	//			}
	//		}
	//	}
	//	END_GRID_ITERATION_SUM(nnz);

	//	BEGIN_GRID_ITERATION_2D(bc.partial_grids[thread_id])
	//	{
	//		if (bc(i, j) > -1)
	//		{
	//			if (bc(i, j) == bc(mid_i, mid_j))
	//			{
	//				nnz_spe++;
	//			}
	//							
	//			// For the boundary update
	//			if (i == i_start)
	//			{
	//				if (bc(i + 1, j) > -1)
	//				{
	//					if (bc(i + 1, j) == bc(mid_i, mid_j))
	//					{
	//						nnz_spe++;
	//					}
	//				}
	//				if (bc(i - 1, j) > -1)
	//				{
	//					if (bc(i - 1, j) == bc(mid_i, mid_j))
	//					{
	//						nnz_spe++;
	//					}
	//				}
	//				if (bc(i, j + 1) > -1)
	//				{
	//					if (bc(i, j + 1) == bc(mid_i, mid_j))
	//					{
	//						nnz_spe++;
	//					}
	//				}
	//				if (bc(i, j - 1) > -1)
	//				{
	//					if (bc(i, j - 1) == bc(mid_i, mid_j))
	//					{
	//						nnz_spe++;
	//					}
	//				} 
	//			}
	//			else if (i == i_end)
	//			{
	//				if (bc(i + 1, j) > -1)
	//				{
	//					if (bc(i + 1, j) == bc(mid_i, mid_j))
	//					{
	//						nnz_spe++;
	//					}
	//				}
	//				if (bc(i - 1, j) > -1)
	//				{
	//					if (bc(i - 1, j) == bc(mid_i, mid_j))
	//					{
	//						nnz_spe++;
	//					}
	//				}
	//				if (bc(i, j + 1) > -1)
	//				{
	//					if (bc(i, j + 1) == bc(mid_i, mid_j))
	//					{
	//						nnz_spe++;
	//					}
	//				}
	//				if (bc(i, j - 1) > -1)
	//				{
	//					if (bc(i, j - 1) == bc(mid_i, mid_j))
	//					{
	//						nnz_spe++;
	//					}
	//				} 
	//			}
	//			else if (j == j_start)
	//			{
	//				if (bc(i + 1, j) > -1)
	//				{
	//					if (bc(i + 1, j) == bc(mid_i, mid_j))
	//					{
	//						nnz_spe++;
	//					}
	//				}
	//				if (bc(i - 1, j) > -1)
	//				{
	//					if (bc(i - 1, j) == bc(mid_i, mid_j))
	//					{
	//						nnz_spe++;
	//					}
	//				}
	//				if (bc(i, j + 1) > -1)
	//				{
	//					if (bc(i, j + 1) == bc(mid_i, mid_j))
	//					{
	//						nnz_spe++;
	//					}
	//				}
	//				if (bc(i, j - 1) > -1)
	//				{
	//					if (bc(i, j - 1) == bc(mid_i, mid_j))
	//					{
	//						nnz_spe++;
	//					}
	//				} 
	//			}
	//			else if (j == j_end)
	//			{
	//				if (bc(i + 1, j) > -1)
	//				{
	//					if (bc(i + 1, j) == bc(mid_i, mid_j))
	//					{
	//						nnz_spe++;
	//					}
	//				}
	//				if (bc(i - 1, j) > -1)
	//				{
	//					if (bc(i - 1, j) == bc(mid_i, mid_j))
	//					{
	//						nnz_spe++;
	//					}
	//				}
	//				if (bc(i, j + 1) > -1)
	//				{
	//					if (bc(i, j + 1) == bc(mid_i, mid_j))
	//					{
	//						nnz_spe++;
	//					}
	//				}
	//				if (bc(i, j - 1) > -1)
	//				{
	//					if (bc(i, j - 1) == bc(mid_i, mid_j))
	//					{
	//						nnz_spe++;
	//					}
	//				} 
	//			}
	//			else
	//			{
	//				if (is_regularized)
	//				{
	//					if (bc(i + 1, j) > -1)
	//					{
	//						nnz++;
	//					}
	//					if (bc(i - 1, j) > -1)
	//					{
	//						nnz++;
	//					}
	//					if (bc(i, j + 1) > -1)
	//					{
	//						nnz++;
	//					}
	//					if (bc(i, j - 1) > -1)
	//					{
	//						nnz++;
	//					}
	//					if (bc(i - 1, j - 1) > -1)
	//					{
	//						nnz++;
	//					}
	//					if (bc(i - 1, j + 1) > -1)
	//					{
	//						nnz++;
	//					}
	//					if (bc(i + 1, j - 1) > -1)
	//					{
	//						nnz++;
	//					}
	//					if (bc(i + 1, j + 1) > -1)
	//					{
	//						nnz++;
	//					} 
	//				}
	//				else
	//				{
	//					if (stencil.values[1].array_for_this(i, j) >= stencil.values[0].array_for_this(i, j))
	//					{
	//						if (bc(i + 1, j) > -1)
	//						{
	//							if (bc(i + 1, j) == bc(mid_i, mid_j))
	//							{
	//								nnz_spe++;
	//							}
	//						}
	//						if (bc(i - 1, j) > -1)
	//						{
	//							if (bc(i - 1, j) == bc(mid_i, mid_j))
	//							{
	//								nnz_spe++;
	//							}
	//						}
	//						if (bc(i, j + 1) > -1)
	//						{
	//							if (bc(i, j + 1) == bc(mid_i, mid_j))
	//							{
	//								nnz_spe++;
	//							}
	//						}
	//						if (bc(i, j - 1) > -1)
	//						{
	//							if (bc(i, j - 1) == bc(mid_i, mid_j))
	//							{
	//								nnz_spe++;
	//							}
	//						} 
	//					}
	//					else
	//					{
	//						if (bc(i - 1, j - 1) > -1)
	//						{
	//							if (bc(i - 1, j - 1) == bc(mid_i, mid_j))
	//							{
	//								nnz_spe++;
	//							}
	//						}
	//						if (bc(i - 1, j + 1) > -1)
	//						{
	//							if (bc(i - 1, j + 1) == bc(mid_i, mid_j))
	//							{
	//								nnz_spe++;
	//							}
	//						}
	//						if (bc(i + 1, j - 1) > -1)
	//						{
	//							if (bc(i + 1, j - 1) == bc(mid_i, mid_j))
	//							{
	//								nnz_spe++;
	//							}
	//						}
	//						if (bc(i + 1, j + 1) > -1)
	//						{
	//							if (bc(i + 1, j + 1) == bc(mid_i, mid_j))
	//							{
	//								nnz_spe++;
	//							}
	//						} 
	//					} 
	//				}
	//			}
	//		}
	//	}
	//	END_GRID_ITERATION_SUM(nnz_spe);

	//	nnz = nnz + (N - nnz_spe) - 3*(bc.grid.i_res);
	//	
	//	cout << nnz << endl;
	//	return nnz;
	//}

	void NewtonMethod(CSR_MATRIX<T>& J, VECTOR_ND<T>& x, VECTOR_ND<T>& b_vector, const int& thread_id)
	{
		// For initializing Jacobian
		const int num_all_full_cells = AssignSequentialindicesToFullCells(bc, thread_id);
		
		BEGIN_HEAD_THREAD_WORK
		{
			multithreading->SplitDomainIndex1D(0, num_all_full_cells);
		}
		END_HEAD_THREAD_WORK;

		BEGIN_HEAD_THREAD_WORK
		{
			x.Initialize(num_all_full_cells, true);	
		}
		END_HEAD_THREAD_WORK;
		
		const int N(x.num_dimension), start_ix(multithreading->start_ix_1D[thread_id]), end_ix(multithreading->end_ix_1D[thread_id]);

		BEGIN_HEAD_THREAD_WORK
		{
			num_iteration = 0;
		}
		END_HEAD_THREAD_WORK;
		
		BEGIN_HEAD_THREAD_WORK
		{
			Jb.Initialize(num_all_full_cells, true);
		}
		END_HEAD_THREAD_WORK;
			
		T* xval(x.values);

		// For Gauss-Seidel Test
		/*CSR_MATRIX<T> A_Matrix;
		A_Matrix.Initialize(3, 9);

		A_Matrix.AssignValue(0, 0, 12);
		A_Matrix.AssignValue(0, 1, 3);
		A_Matrix.AssignValue(0, 2, -5);
		A_Matrix.AssignValue(1, 0, 1);
		A_Matrix.AssignValue(1, 1, 5);
		A_Matrix.AssignValue(1, 2, 3);
		A_Matrix.AssignValue(2, 0, 3);
		A_Matrix.AssignValue(2, 1, 7);
		A_Matrix.AssignValue(2, 2, 13);

		VECTOR_ND<T> b_vector_test;
		b_vector_test.Initialize(3, true);
		
		b_vector_test[0] = 1;
		b_vector_test[1] = 28;
		b_vector_test[2] = 76;
		
		VECTOR_ND<T> x_vector_sol;
		x_vector_sol.Initialize(3, true);

		x_vector_sol[0] = 1;
		x_vector_sol[1] = 0;
		x_vector_sol[2] = 1;

		sub_linear_solver->Solve(A_Matrix, x_vector_sol, b_vector_test, bc, thread_id);

		cout << x_vector_sol[0] << endl;
		cout << x_vector_sol[1] << endl;
		cout << x_vector_sol[2] << endl;*/

		T stopping_criterion(0);
		T temp;

		// For Newton Method Test
		//VECTOR_ND<T> newton_test;
		//newton_test.Initialize(3, true);
		//newton_test[0] = 1;
		//newton_test[1] = 0;
		//newton_test[2] = 1;
		//CSR_MATRIX<T> test_jacobian;
		//test_jacobian.Initialize(3, 9);
		//VECTOR_ND<T> test_JB;
		//test_JB.Initialize(3, true);
		//VECTOR_3D<T> test_b_vector;
		//		
		//while (num_iteration < max_iteration_for_Newton_Method)
		//{
		//	test_b_vector.x = POW2(newton_test[0]) + POW2(newton_test[1]) + POW2(newton_test[2]) - 3;
		//	test_b_vector.y = POW2(newton_test[0]) + POW2(newton_test[1]) - newton_test[2] - 1;
		//	test_b_vector.z = newton_test[0] + newton_test[1] + newton_test[2] - 3;

		//	cout << test_b_vector.x << endl;
		//	cout << test_b_vector.y << endl;
		//	cout << test_b_vector.z << endl;

		//	test_jacobian.AssignValue(0, 0, 2*newton_test[0]);
		//	test_jacobian.AssignValue(0, 1, 2*newton_test[1]);
		//	test_jacobian.AssignValue(0, 2, 2*newton_test[2]);
		//	test_jacobian.AssignValue(1, 0, 2*newton_test[0]);
		//	test_jacobian.AssignValue(1, 1, 2*newton_test[1]);
		//	test_jacobian.AssignValue(1, 2, -1);
		//	test_jacobian.AssignValue(2, 0, 1);
		//	test_jacobian.AssignValue(2, 1, 1);
		//	test_jacobian.AssignValue(2, 2, 1);
		//	
		//	MATRIX_3X3 for_inverse(2*newton_test[0], 2*newton_test[1], 2*newton_test[2], 2*newton_test[0], 2*newton_test[1], -1, 1, 1, 1);

		//	BEGIN_HEAD_THREAD_WORK
		//	{
		//		temp = stopping_criterion;
		//	}
		//	END_HEAD_THREAD_WORK;
		//				
		//	//sub_linear_solver->Solve(test_jacobian, test_JB, test_b_vector, bc, thread_id);

		//	Jb[0] = (for_inverse.Inversed()*test_b_vector)[0];
		//	Jb[1] = (for_inverse.Inversed()*test_b_vector)[1];
		//	Jb[2] = (for_inverse.Inversed()*test_b_vector)[2];
		//	
		//	for (int i = 0; i <= 2; i++)
		//	{
		//		newton_test[i] = newton_test[i] - alpha*Jb[i];
		//	}
		//	multithreading->Sync(thread_id);
		//	
		//	// absolute max-norm residual
		//	/*T max_norm(0);
		//	for (int i = start_ix; i <= end_ix; i++)
		//	{
		//		max_norm = MAX(max_norm,abs(alpha*Jb[i]));
		//	}
		//	multithreading->SyncMax(thread_id, max_norm);
		//	
		//	stopping_criterion = max_norm;*/

		//	// Discrete 1-norm
		//	for (int i = 0; i <= 2; i++)
		//	{
		//		stopping_criterion += abs(alpha*Jb[i]);
		//	}
		//	multithreading->SyncSum(thread_id, stopping_criterion);

		//	BEGIN_HEAD_THREAD_WORK
		//	{
		//		num_iteration++;
		//		cout << "Increment: " << stopping_criterion << endl;
		//		cout << "alpha: " << alpha << endl;
		//		cout << newton_test[0] << endl;
		//		cout << newton_test[1] << endl;
		//		cout << newton_test[2] << endl;
		//	}
		//	END_HEAD_THREAD_WORK;

		//	// absolute 2-norm residual
		//	/*for (int i = start_ix; i <= end_ix; i++)
		//	{
		//		stopping_criterion += POW2(alpha*Jb[i]);
		//	}
		//	multithreading->SyncSum(thread_id, stopping_criterion);

		//	stopping_criterion = sqrt(stopping_criterion);*/

		//	//VectorToGrid(Jb, u, bc, thread_id);
		//	/*ofstream fout;
		//	fout.open("Jb");
		//	for (int j = u.j_start; j <= u.j_end; j++)
		//	{
		//		for (int i = u.i_start; i <= u.i_end; i++)
		//			fout << u(i, j) << " ";

		//			fout << "\n";
		//	}
		//	fout.close();*/

		//	//if (stopping_criterion < tolerance)
		//	//{
		//	//	cout << "--------------Newton Method--------------" << endl;
		//	//	cout << "Converge!!" << endl;
		//	//	cout << "Iteration Number : " << num_iteration << endl;
		//	//	cout << "Increment: " << stopping_criterion << endl;
		//	//	cout << "-----------------------------------------" << endl;
		//	//	
		//	//	/*ofstream fout;
		//	//	fout.open("solution");
		//	//	for (int j = u.j_start; j <= u.j_end; j++)
		//	//	{
		//	//		for (int i = u.i_start; i <= u.i_end; i++)
		//	//			fout << u(i, j) << " ";

		//	//			fout << "\n";
		//	//	}
		//	//	fout.close();*/
		//	//	break;
		//	//}
		//	//if (num_iteration == max_iteration_for_Newton_Method - 1)
		//	//{
		//	//	cout << "--------------Newton Method--------------" << endl;
		//	//	cout << "Not Converge!!" << endl;
		//	//	cout << "Iteration Number : " << num_iteration << endl;
		//	//	cout << "Increment: " << stopping_criterion << endl;
		//	//	cout << "-----------------------------------------" << endl;
		//	//}

		//	//BEGIN_HEAD_THREAD_WORK
		//	//{
		//	//	num_iteration++;
		//	//	cout << "Increment: " << stopping_criterion << endl;
		//	//	cout << "alpha: " << alpha << endl;
		//	//	
		//	//	if (temp <= stopping_criterion)
		//	//	{
		//	//		alpha = (T)0.5*alpha;
		//	//	}
		//	//}
		//	//END_HEAD_THREAD_WORK;
		//}

		while (num_iteration < max_iteration_for_Newton_Method)
		{
			BEGIN_HEAD_THREAD_WORK
			{
				temp = stopping_criterion;
			}
			END_HEAD_THREAD_WORK;

			BEGIN_HEAD_THREAD_WORK
			{
				stopping_criterion = 0;
			}
			END_HEAD_THREAD_WORK;
			
			SetupDensity(thread_id);
			
			CalculateJacobian(J, b_vector, thread_id);
			
			/*ofstream fout;
			fout.open("one_over_2dx");
			fout << bc.grid.one_over_2dx << endl;
			fout.close();*/

			/*fout.open("val");
			for (int i = 0; i < J.nz; i++)
			{
				fout << J.values[i] << endl;
			}
			fout.close();

			fout.open("ci");
			for (int i = 0; i < J.nz; i++)
			{
				fout << J.column_index[i] << endl;
			}
			fout.close();

			fout.open("rp");
			for (int i = 0; i <= J.N; i++)
			{
				fout << J.row_ptr[i] << endl;
			}
			fout.close();

			fout.open("b_vector");
			for (int i = 0; i < J.N; i++)
			{
				fout << b_vector[i] << endl;
			}
			fout.close();*/

			sub_linear_solver->Solve(J, Jb, b_vector, bc, thread_id);

			/*fout.open("Jb_vector");
			for (int i = 0; i < J.N; i++)
			{
				fout << Jb[i] << endl;
			}
			fout.close();*/

			GridToVector(u, x, bc, thread_id);

			/*fout.open("u_vector");
			for (int i = 0; i < J.N; i++)
			{
				fout << x[i] << endl;
			}
			fout.close();

			fout.open("bc");
			for (int i = bc.i_start; i <= bc.i_end; i++)
			{
				for (int j = bc.j_start; j <= bc.j_end; j++)
				{
					fout << bc(i, j) << " ";
				}
				fout << endl;
			}
			fout.close();*/

			for (int i = start_ix; i <= end_ix; i++)
			{
				xval[i] = xval[i] - alpha*Jb[i];
				
			}
			multithreading->Sync(thread_id);
			
			// absolute max-norm residual
			/*T max_norm(0);
			for (int i = start_ix; i <= end_ix; i++)
			{
				max_norm = MAX(max_norm,abs(alpha*Jb[i]));
			}
			multithreading->SyncMax(thread_id, max_norm);
			
			stopping_criterion = max_norm;*/

			// Discrete 1-norm
			for (int i = start_ix; i <= end_ix; i++)
			{
				stopping_criterion += u.grid.dx*u.grid.dx*abs(alpha*Jb[i]);
			}
			multithreading->SyncSum(thread_id, stopping_criterion);

			// absolute 2-norm residual
			/*for (int i = start_ix; i <= end_ix; i++)
			{
				stopping_criterion += POW2(alpha*Jb[i]);
			}
			multithreading->SyncSum(thread_id, stopping_criterion);

			stopping_criterion = sqrt(stopping_criterion);*/

			VectorToGrid(x, u, bc, thread_id);

			//VectorToGrid(Jb, u, bc, thread_id);
			/*ofstream fout;
			fout.open("Jb");
			for (int j = u.j_start; j <= u.j_end; j++)
			{
				for (int i = u.i_start; i <= u.i_end; i++)
					fout << u(i, j) << " ";

					fout << "\n";
			}
			fout.close();*/

			if (stopping_criterion < tolerance)
			{
				cout << "--------------Newton Method--------------" << endl;
				cout << "Converge!!" << endl;
				cout << "Iteration Number : " << num_iteration << endl;
				cout << "Increment: " << stopping_criterion << endl;
				cout << "-----------------------------------------" << endl;
				
				/*ofstream fout;
				fout.open("solution");
				for (int j = u.j_start; j <= u.j_end; j++)
				{
					for (int i = u.i_start; i <= u.i_end; i++)
						fout << u(i, j) << " ";

						fout << "\n";
				}
				fout.close();

				fout.open("true_solution_x");
				for (int j = true_solution.j_start; j <= true_solution.j_end; j++)
				{
					for (int i = true_solution.i_start; i <= true_solution.i_end; i++)
					{
						fout << true_solution(i, j).x << " ";
					}
					fout << "\n";
				}
				fout.close();

				fout.open("true_solution_y");
				for (int j = true_solution.j_start; j <= true_solution.j_end; j++)
				{
					for (int i = true_solution.i_start; i <= true_solution.i_end; i++)
					{
						fout << true_solution(i, j).y << " ";
					}
					fout << "\n";
				}
				fout.close();*/

				break;
			}
			if (num_iteration == max_iteration_for_Newton_Method - 1)
			{
				cout << "--------------Newton Method--------------" << endl;
				cout << "Not Converge!!" << endl;
				cout << "Iteration Number : " << num_iteration << endl;
				cout << "Increment: " << stopping_criterion << endl;
				cout << "-----------------------------------------" << endl;
				
				/*ofstream fout;
				fout.open("solution");
				for (int j = u.j_start; j <= u.j_end; j++)
				{
					for (int i = u.i_start; i <= u.i_end; i++)
						fout << u(i, j) << " ";

						fout << "\n";
				}
				fout.close();

				fout.open("true_solution");
				for (int j = true_solution.j_start; j <= true_solution.j_end; j++)
				{
					for (int i = true_solution.i_start; i <= true_solution.i_end; i++)
					{
						fout << true_solution(i, j) << " ";
					}
					fout << "\n";
				}*/
			}

			BEGIN_HEAD_THREAD_WORK
			{
				num_iteration++;
				cout << "Increment: " << stopping_criterion << endl;
				cout << "alpha: " << alpha << endl;

				if (temp <= stopping_criterion)
				{
					alpha = (T)0.5*alpha;
				}
			}
			END_HEAD_THREAD_WORK;
		}
	}

	int AssignSequentialindicesToFullCells(const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id)
	{
		// Count number of full cells
		int full_ix(0);

		BEGIN_GRID_ITERATION_2D(bc.partial_grids[thread_id])
		{
			if (bc(i, j) > -1)
			{
				++full_ix;
			}	
		}
		END_GRID_ITERATION_2D;

		// Indexing the value of boundary condition field
		int start_full_ix, end_full_ix;
		multithreading->SyncDomainIndices1D(thread_id, full_ix, start_full_ix, end_full_ix);

		BEGIN_GRID_ITERATION_2D(bc.partial_grids[thread_id])
		{
			if (bc(i, j) > -1)
			{
				bc(i, j) = start_full_ix++;
			}		
		}
		END_GRID_ITERATION_2D;

		assert(start_full_ix - 1 == end_full_ix);

		multithreading->SyncSum(thread_id, full_ix);

		return full_ix;
	}

	void SetupBoundaryCondition(FIELD_STRUCTURE_2D<int>& bc_input, const int& thread_id)
	{
		ARRAY_2D<int>& bc_array(bc_input.array_for_this);
		GRID_STRUCTURE_2D& grid(bc_input.grid);

		int i(0), j(0);
		
		BEGIN_GRID_ITERATION_2D(bc_input.partial_grids_ghost[thread_id])
		{
			// Speed-up variable
			if (i < grid.i_start || i > grid.i_end || j < grid.j_start || j > grid.j_end)
			{
				bc_array(i, j) = BC_IMPLICIT;
			}
			else
			{
				bc_array(i, j) = BC_FULL;
			}
		}
		END_GRID_ITERATION_2D;
	}

	void SetupDensity(const int& thread_id)
	{
		if (test_number == 1)
		{
			ARRAY<T> sub(density_x.grid.i_res);
			ARRAY<T> sub_x(density_x.grid.i_res);
			ARRAY<T> sub_xx(density_x.grid.i_res);
			
			// For the integration
			ARRAY<T> sub_h(density_x.grid.i_res);
			ARRAY<T> sub_h_x(density_x.grid.i_res);
			ARRAY<T> sub_h_xx(density_x.grid.i_res);

			BEGIN_HEAD_THREAD_WORK
			{
				multithreading->SplitDomainIndex1D(0, sub.length);
			}
			END_HEAD_THREAD_WORK

			const int N(sub.length), start_ix(multithreading->start_ix_1D[thread_id]), end_ix(multithreading->end_ix_1D[thread_id]);

			for (int i = start_ix; i <= end_ix; i++)
			{
				T coor = density_x.x_min + i*density_x.dx;
				sub[i] = (-(T)1/(8*PI)*POW2(coor) + (T)1/(256*PI*PI*PI) + (T)1/(32*PI))*cos(8*PI*coor) + (T)1/(32*PI*PI)*coor*sin(8*PI*coor);
				sub_x[i] = -(T)2/(8*PI)*coor*cos(8*PI*coor) -8*PI*(-(T)1/(8*PI)*POW2(coor) + (T)1/(256*PI*PI*PI) + (T)1/(32*PI))*sin(8*PI*coor) + (T)1/(32*PI*PI)*sin(8*PI*coor) + (T)8*PI/(32*PI*PI)*coor*cos(8*PI*coor);
				sub_xx[i] = -(T)2/(8*PI)*cos(8*PI*coor) + (T)2/(8*PI)*8*PI*coor*sin(8*PI*coor) - 64*PI*PI*(-(T)1/(8*PI)*POW2(coor) + (T)1/(256*PI*PI*PI) + (T)1/(32*PI))*cos(8*PI*coor) - (8*PI)*(-(T)2/(8*PI)*coor)*sin(8*PI*coor) + (T)8*PI/(32*PI*PI)*cos(8*PI*coor) + (T)8*PI/(32*PI*PI)*cos(8*PI*coor) - (T)64*PI*PI/(32*PI*PI)*coor*sin(8*PI*coor);
			}
			multithreading->Sync(thread_id);

			BEGIN_GRID_ITERATION_2D(density_x.partial_grids[thread_id])
			{
				density_x(i, j) = 1 + 4*(sub_xx[i]*sub[j] + sub[i]*sub_xx[j]) + 16*(sub[i]*sub[j]*sub_xx[i]*sub_xx[j] - POW2(sub_x[i])*POW2(sub_x[j]));
			}
			END_GRID_ITERATION_2D;

			BEGIN_GRID_ITERATION_2D(density_y.partial_grids[thread_id])
			{
				density_y(i, j) = (T)1;
			}
			END_GRID_ITERATION_2D;
			
			// Define partial_p, partial_q
			BEGIN_GRID_ITERATION_2D(partial_p.partial_grids[thread_id])
			{
				partial_p(i, j) = 0;
			}
			END_GRID_ITERATION_2D;

			BEGIN_GRID_ITERATION_2D(partial_q.partial_grids[thread_id])
			{
				partial_q(i, j) = 0;
			}
			END_GRID_ITERATION_2D;

			// True solution
			BEGIN_GRID_ITERATION_2D(true_solution.partial_grids[thread_id])
			{
				true_solution(i, j).x = (true_solution.grid.x_min + i*true_solution.grid.dx) + 4*sub_x[i]*sub[j];
				true_solution(i, j).y = (true_solution.grid.y_min + j*true_solution.grid.dy) + 4*sub[i]*sub_x[j];
			}
			END_GRID_ITERATION_2D;
		}
		
		if (test_number == 2)
		{
			MATRIX_2X2 M_x(0.8, 0, 0, 0.4), M_y(0.6, 0.2, 0.2, 0.8), J(0, 1, -1, 0);

			MATRIX_2X2 nu, deno;

			nu = M_x.Inversed()*M_y.Inversed()*J;
			deno = M_x.Inversed()*M_y.Inversed();
			
			T atheta = nu.Trace()/deno.Trace();

			MATRIX_2X2 R_theta(cos(atan(atheta)), sin(atan(atheta)), -sin(atan(atheta)), cos(atan(atheta)));

			MATRIX_2X2 Ellipse_Mapping;

			Ellipse_Mapping = M_y*R_theta*M_x.Inversed();

			BEGIN_GRID_ITERATION_2D(density_x.partial_grids[thread_id])
			{
				//density_x(i, j) = Ellipse_Mapping.x[0]*Ellipse_Mapping.x[3] - POW2(Ellipse_Mapping.x[2]);
				T x_coor = density_x.grid.x_min + i*density_x.grid.dx, y_coor = density_x.grid.y_min + j*density_y.grid.dy;
				T domain_region = POW2(x_coor)/POW2(0.8) + POW2(y_coor)/POW2(0.4);
				
				if (domain_region < 1)
				{
					density_x(i, j) = Ellipse_Mapping.x[0]*Ellipse_Mapping.x[3] - POW2(Ellipse_Mapping.x[2]);
				}
				else
				{
					density_x(i, j) = (T)0;
				}
			}
			END_GRID_ITERATION_2D;

			BEGIN_GRID_ITERATION_2D(density_y.partial_grids[thread_id])
			{
				density_y(i, j) = (T)1;
			}
			END_GRID_ITERATION_2D;

			// True solution
			BEGIN_GRID_ITERATION_2D(true_solution.partial_grids[thread_id])
			{
				T x_coor = true_solution.grid.x_min + i*true_solution.grid.dx, y_coor = true_solution.grid.y_min + j*true_solution.grid.dy;
				
				true_solution(i, j).x = Ellipse_Mapping.x[0]*x_coor + Ellipse_Mapping.x[1]*y_coor;
				true_solution(i, j).y = Ellipse_Mapping.x[2]*x_coor + Ellipse_Mapping.x[3]*y_coor;

				if (POW2(x_coor)/POW2(0.8) + POW2(y_coor)/POW2(0.4) <= 1)
				{
					true_solution.draw_for_this(i, j) = true;
				}
				else
				{
					true_solution.draw_for_this(i, j) = false;
				}

				if (POW2(x_coor)/POW2(0.8) + POW2(y_coor)/POW2(0.4) <= 1)
				{
					grad_u.draw_for_this(i, j) = true;
				}
				else
				{
					grad_u.draw_for_this(i, j) = false;
				}
			}
			END_GRID_ITERATION_2D;
		}

		if (test_number == 3)
		{
			BEGIN_GRID_ITERATION_2D(density_x.partial_grids[thread_id])
			{
				T x_coor = density_x.grid.x_min + i*density_x.grid.dx, y_coor = density_x.grid.y_min + j*density_y.grid.dy;
				
				T domain_region_1 = POW2(x_coor + 0.1) + POW2(y_coor);
				T domain_region_2 = POW2(x_coor - 0.1) + POW2(y_coor);
				
				if (domain_region_1 < POW2(0.85) && x_coor < -0.1)
				{
					density_x(i, j) = 1;
				}
				else if (domain_region_2 < POW2(0.85) && x_coor > 0.1)
				{
					density_x(i, j) = 1;
				}
				else
				{
					density_x(i, j) = (T)0;
				}
			}
			END_GRID_ITERATION_2D;

			BEGIN_GRID_ITERATION_2D(density_y.partial_grids[thread_id])
			{
				density_y(i, j) = (T)1;
			}
			END_GRID_ITERATION_2D;

			// True solution
			BEGIN_GRID_ITERATION_2D(true_solution.partial_grids[thread_id])
			{
				T x_coor = true_solution.grid.x_min + i*true_solution.grid.dx, y_coor = true_solution.grid.y_min + j*true_solution.grid.dy;
				
				T domain_region_1 = POW2(x_coor + 0.1) + POW2(y_coor);
				T domain_region_2 = POW2(x_coor - 0.1) + POW2(y_coor);

				if (x_coor < -0.1)
				{
					true_solution(i, j).x = x_coor + 0.1;
					true_solution(i, j).y = y_coor;	
				}
				if (x_coor > 0.1)
				{
					true_solution(i, j).x = x_coor - 0.1;
					true_solution(i, j).y = y_coor;	
				}

				if (domain_region_1 < POW2(0.85) && x_coor < -0.1)
				{
					true_solution.draw_for_this(i, j) = true;
				}
				else if (domain_region_2 < POW2(0.85) && x_coor > 0.1)
				{
					true_solution.draw_for_this(i, j) = true;
				}
				else
				{
					true_solution.draw_for_this(i, j) = false;
				}

				if (domain_region_1 <= POW2(0.85) && x_coor <= -0.1)
				{
					grad_u.draw_for_this(i, j) = true;
				}
				else if (domain_region_2 <= POW2(0.85) && x_coor >= 0.1)
				{
					grad_u.draw_for_this(i, j) = true;
				}
				else
				{
					grad_u.draw_for_this(i, j) = false;
				}
			}
			END_GRID_ITERATION_2D;
		}
		
		// Prins's thesis - Example 7.4.1
		if (test_number == 8)
		{
			T deno = (T)4*0.247994810805468;

			T G_0 = (T)1/deno;

			/*T deno = POW2(PI/8) + 1;

			T G_0 = ((T)1/PI/4)/(1 - (T)1/deno);*/

			/*T a = cos(5*PI/8), b = cos(3*PI/8);

			T deno_1 = POW2((1 + b)/(1 - b)) + 1, deno_2 = POW2((1 + a)/(1 - a)) + 1;

			T G_0 = (T)1/(4*PI*(-(T)1/deno_1 + (T)1/deno_2));*/

			//T a = cos(PI/8);

			//T deno = (1 + a)/(1 - a) + 1;

			//T G_0 = 10;//(T)1/(4*PI)*deno;

			/*T a = sin(PI/8);

			T b = sqrt(POW2((T)1/a) - 1) + (T)1/a;

			T deno = 2*PI*(POW2(b) + 1);

			T G_0 = (T)1/deno;*/

			BEGIN_GRID_ITERATION_2D(density_x.partial_grids[thread_id])
			{
				density_x(i, j) = (T)1;
			}
			END_GRID_ITERATION_2D;

			BEGIN_GRID_ITERATION_2D(density_y.partial_grids[thread_id])
			{
				T p_coor, q_coor;
				
				if (i == i_start)
				{
					if (j == j_start)
					{
						p_coor = (u(i + 1, j) - u(i, j))*density_y.one_over_dx, q_coor = (u(i, j + 1) - u(i, j))*density_y.one_over_dy;
					}
					else if (j == j_end)
					{
						p_coor = (u(i + 1, j) - u(i, j))*density_y.one_over_dx, q_coor = (u(i, j) - u(i, j - 1))*density_y.one_over_dy;
					}
					else
					{
						p_coor = (u(i + 1, j) - u(i, j))*density_y.one_over_dx, q_coor = (u(i, j + 1) - u(i, j - 1))*density_y.one_over_2dy;
					}
				}
				else if (i == i_end)
				{
					if (j == j_start)
					{
						p_coor = (u(i, j) - u(i - 1, j))*density_y.one_over_dx, q_coor = (u(i, j + 1) - u(i, j))*density_y.one_over_dy;
					}
					else if (j == j_end)
					{
						p_coor = (u(i, j) - u(i - 1, j))*density_y.one_over_dx, q_coor = (u(i, j) - u(i, j - 1))*density_y.one_over_dy;
					}
					else
					{
						p_coor = (u(i, j) - u(i - 1, j))*density_y.one_over_dx, q_coor = (u(i, j + 1) - u(i, j - 1))*density_y.one_over_2dy;
					}
				}
				else if (j == j_start)
				{
					p_coor = (u(i + 1, j) - u(i - 1, j))*density_y.one_over_2dx, q_coor = (u(i, j + 1) - u(i, j))*density_y.one_over_dy;
				}
				else if (j == j_end)
				{
					p_coor = (u(i + 1, j) - u(i - 1, j))*density_y.one_over_2dx, q_coor = (u(i, j) - u(i, j - 1))*density_y.one_over_dy;
				}
				else
				{
					p_coor = (u(i + 1, j) - u(i - 1, j))*density_y.one_over_2dx, q_coor = (u(i, j + 1) - u(i, j - 1))*density_y.one_over_2dy;
				}
				
				density_y(i, j) = 4*G_0/POW2(POW2(p_coor) + POW2(q_coor) + 1);
			}
			END_GRID_ITERATION_2D;

			// Define partial_p, partial_q
			BEGIN_GRID_ITERATION_2D(partial_p.partial_grids[thread_id])
			{
				T p_coor, q_coor;
				
				if (i == i_start)
				{
					if (j == j_start)
					{
						p_coor = (u(i + 1, j) - u(i, j))*density_y.one_over_dx, q_coor = (u(i, j + 1) - u(i, j))*density_y.one_over_dy;
					}
					else if (j == j_end)
					{
						p_coor = (u(i + 1, j) - u(i, j))*density_y.one_over_dx, q_coor = (u(i, j) - u(i, j - 1))*density_y.one_over_dy;
					}
					else
					{
						p_coor = (u(i + 1, j) - u(i, j))*density_y.one_over_dx, q_coor = (u(i, j + 1) - u(i, j - 1))*density_y.one_over_2dy;
					}
					
				}
				else if (i == i_end)
				{
					if (j == j_start)
					{
						p_coor = (u(i, j) - u(i - 1, j))*density_y.one_over_dx, q_coor = (u(i, j + 1) - u(i, j))*density_y.one_over_dy;
					}
					else if (j == j_end)
					{
						p_coor = (u(i, j) - u(i - 1, j))*density_y.one_over_dx, q_coor = (u(i, j) - u(i, j - 1))*density_y.one_over_dy;
					}
					else
					{
						p_coor = (u(i, j) - u(i - 1, j))*density_y.one_over_dx, q_coor = (u(i, j + 1) - u(i, j - 1))*density_y.one_over_2dy;
					}
				}
				else if (j == j_start)
				{
					p_coor = (u(i + 1, j) - u(i - 1, j))*density_y.one_over_2dx, q_coor = (u(i, j + 1) - u(i, j))*density_y.one_over_dy;
				}
				else if (j == j_end)
				{
					p_coor = (u(i + 1, j) - u(i - 1, j))*density_y.one_over_2dx, q_coor = (u(i, j) - u(i, j - 1))*density_y.one_over_dy;
				}
				else
				{
					p_coor = (u(i + 1, j) - u(i - 1, j))*density_y.one_over_2dx, q_coor = (u(i, j + 1) - u(i, j - 1))*density_y.one_over_2dy;
				}

				partial_p(i, j) = (T)1/G_0*(POW2(p_coor) + POW2(q_coor) + 1)*p_coor;
			}
			END_GRID_ITERATION_2D;

			BEGIN_GRID_ITERATION_2D(partial_q.partial_grids[thread_id])
			{
				T p_coor, q_coor;
				
				if (i == i_start)
				{
					if (j == j_start)
					{
						p_coor = (u(i + 1, j) - u(i, j))*density_y.one_over_dx, q_coor = (u(i, j + 1) - u(i, j))*density_y.one_over_dy;
					}
					else if (j == j_end)
					{
						p_coor = (u(i + 1, j) - u(i, j))*density_y.one_over_dx, q_coor = (u(i, j) - u(i, j - 1))*density_y.one_over_dy;
					}
					else
					{
						p_coor = (u(i + 1, j) - u(i, j))*density_y.one_over_dx, q_coor = (u(i, j + 1) - u(i, j - 1))*density_y.one_over_2dy;
					}
					
				}
				else if (i == i_end)
				{
					if (j == j_start)
					{
						p_coor = (u(i, j) - u(i - 1, j))*density_y.one_over_dx, q_coor = (u(i, j + 1) - u(i, j))*density_y.one_over_dy;
					}
					else if (j == j_end)
					{
						p_coor = (u(i, j) - u(i - 1, j))*density_y.one_over_dx, q_coor = (u(i, j) - u(i, j - 1))*density_y.one_over_dy;
					}
					else
					{
						p_coor = (u(i, j) - u(i - 1, j))*density_y.one_over_dx, q_coor = (u(i, j + 1) - u(i, j - 1))*density_y.one_over_2dy;
					}
				}
				else if (j == j_start)
				{
					p_coor = (u(i + 1, j) - u(i - 1, j))*density_y.one_over_2dx, q_coor = (u(i, j + 1) - u(i, j))*density_y.one_over_dy;
				}
				else if (j == j_end)
				{
					p_coor = (u(i + 1, j) - u(i - 1, j))*density_y.one_over_2dx, q_coor = (u(i, j) - u(i, j - 1))*density_y.one_over_dy;
				}
				else
				{
					p_coor = (u(i + 1, j) - u(i - 1, j))*density_y.one_over_2dx, q_coor = (u(i, j + 1) - u(i, j - 1))*density_y.one_over_2dy;
				}

				partial_q(i, j) = (T)1/G_0*(POW2(p_coor) + POW2(q_coor) + 1)*q_coor;
			}
			END_GRID_ITERATION_2D;
		}
	}

	void SetupInitialForNewtonMethod(const int& thread_id)
	{
		T a_min, a_max, b_min, b_max, c_min, c_max, d_min, d_max;

		if (test_number == 1)
		{
			a_min = u.grid.x_min, a_max = u.grid.x_max;
			b_min = u.grid.y_min, b_max = u.grid.y_max;
			c_min = u.grid.x_min, c_max = u.grid.x_max;
			d_min = u.grid.y_min, d_max = u.grid.y_max;

			BEGIN_GRID_ITERATION_2D(u.partial_grids[thread_id])
			{
				T x_coor = u.grid.x_min + i*u.grid.dx, y_coor = u.grid.y_min + j*u.grid.dy;

				u(i, j) = (T)0.5*(c_max - c_min)/(a_max - a_min)*POW2(x_coor) + (c_min*a_max - c_max*a_min)/(a_max - a_min)*x_coor  + (T)0.5*(d_max - d_min)/(b_max - b_min)*POW2(y_coor) + (d_min*b_max - d_max*b_min)/(b_max - b_min)*y_coor;
			}
			END_GRID_ITERATION_2D;
			
			//SetupBoundaryCondition(bc, thread_id);

			//BEGIN_GRID_ITERATION_2D(RHS.partial_grids[thread_id])
			//{
			//	RHS(i, j) = sqrt(density_x(i, j)/density_y(i, j));
			//}
			//END_GRID_ITERATION_2D;

			//// For integration
			//BEGIN_GRID_ITERATION_2D(RHS_h.partial_grids[thread_id])
			//{
			//	RHS_h(i, j) = sqrt(density_h_x(i, j)/density_h_y(i, j));
			//}
			//END_GRID_ITERATION_2D;

			//T integral_value(0);

			//BEGIN_GRID_ITERATION_2D(RHS_h.partial_grids[thread_id])
			//{
			//	integral_value += RHS_h(i, j)*RHS_h.grid.dx*RHS_h.grid.dy;
			//}
			//END_GRID_ITERATION_SUM(integral_value);	
			//
			//T rhs_value(0);

			//BEGIN_GRID_ITERATION_2D(RHS_h.partial_grids[thread_id])
			//{
			//	rhs_value += 2*RHS_h.grid.dx*RHS_h.grid.dy;
			//}
			//END_GRID_ITERATION_SUM(rhs_value);

			//T k = rhs_value/integral_value;

			//BEGIN_GRID_ITERATION_2D(RHS.partial_grids[thread_id])
			//{
			//	RHS(i, j) = k*RHS(i, j);
			//}
			//END_GRID_ITERATION_2D;
			//
			//poisson_solver.Solve(u, bc, RHS, thread_id);
		}
		if (test_number == 2)
		{
			BEGIN_GRID_ITERATION_2D(u.partial_grids[thread_id])
			{
				T x_coor = u.grid.x_min + i*u.grid.dx, y_coor = u.grid.y_min + j*u.grid.dy;
				T theta_1 = atan(-(T)1/3), theta_2 = atan(-4);

				a_min = u.grid.x_min, a_max = u.grid.x_max;
				b_min = u.grid.y_min, b_max = u.grid.y_max;
				c_min = -(0.6*cos(theta_1) - 0.2*sin(theta_1)), c_max = -c_min;
				d_min = -(0.2*cos(theta_2) - 0.8*sin(theta_2)), d_max = -d_min;
				
				u(i, j) = (T)0.5*(c_max - c_min)/(a_max - a_min)*POW2(x_coor) + (c_min*a_max - c_max*a_min)/(a_max - a_min)*x_coor  + (T)0.5*(d_max - d_min)/(b_max - b_min)*POW2(y_coor) + (d_min*b_max - d_max*b_min)/(b_max - b_min)*y_coor;
			}
			END_GRID_ITERATION_2D;
		}
		if (test_number == 3)
		{
			BEGIN_GRID_ITERATION_2D(u.partial_grids[thread_id])
			{
				T x_coor = u.grid.x_min + i*u.grid.dx, y_coor = u.grid.y_min + j*u.grid.dy;

				a_min = u.grid.x_min, a_max = u.grid.x_max;
				b_min = u.grid.y_min, b_max = u.grid.y_max;
				c_min = -0.85, c_max = -c_min;
				d_min = -0.85, d_max = -d_min;
				
				u(i, j) = (T)0.5*(c_max - c_min)/(a_max - a_min)*POW2(x_coor) + (c_min*a_max - c_max*a_min)/(a_max - a_min)*x_coor  + (T)0.5*(d_max - d_min)/(b_max - b_min)*POW2(y_coor) + (d_min*b_max - d_max*b_min)/(b_max - b_min)*y_coor;
			}
			END_GRID_ITERATION_2D;
		}
		if (test_number == 8)
		{
			BEGIN_GRID_ITERATION_2D(u.partial_grids[thread_id])
			{
				T x_coor = u.grid.x_min + i*u.grid.dx, y_coor = u.grid.y_min + j*u.grid.dy;

				T deno = cos(PI/8);

				a_min = u.grid.x_min, a_max = u.grid.x_max;
				b_min = u.grid.y_min, b_max = u.grid.y_max;
				c_min = -sqrt((T)1/POW2(deno) - 1), c_max = -c_min;
				d_min = (T)1/deno - sqrt((T)1/POW2(deno) - 1), d_max = (T)1/deno + sqrt((T)1/POW2(deno) - 1);
				
				u(i, j) = (T)0.5*(c_max - c_min)/(a_max - a_min)*POW2(x_coor) + (c_min*a_max - c_max*a_min)/(a_max - a_min)*x_coor  + (T)0.5*(d_max - d_min)/(b_max - b_min)*POW2(y_coor) + (d_min*b_max - d_max*b_min)/(b_max - b_min)*y_coor;
			}
			END_GRID_ITERATION_2D;
		}
	}

	void VectorToGrid(const VECTOR_ND<T>& x, FIELD_STRUCTURE_2D<T>& solution, const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id)
	{
		BEGIN_GRID_ITERATION_2D(solution.partial_grids[thread_id])
		{
			solution.array_for_this(i, j) = x[bc(i, j)];
		}
		END_GRID_ITERATION_2D;
	}

	void GridToVector(const FIELD_STRUCTURE_2D<T>& solution, VECTOR_ND<T>& x, const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id)
	{
		BEGIN_GRID_ITERATION_2D(solution.partial_grids[thread_id])
		{
			x[bc(i, j)] = solution.array_for_this(i, j);
		}
		END_GRID_ITERATION_2D;
	}

	void ComputeGradient(void)
	{
		// Speedup variables
		int i_start(u.grid.i_start), i_end(u.grid.i_end), j_start(u.grid.j_start), j_end(u.grid.j_end);

		for (int i = i_start; i <= i_end; i++)
		{
			for (int j = j_start; j <= j_end; j++)
			{
				VT &grad(grad_u(i, j));

				T grad_x, grad_y;
				
				if (i == i_start)
				{
					if (j == j_start)
					{
						grad_x = (u(i + 1, j) - u(i, j))*u.grid.one_over_dx, grad_y = (u(i, j + 1) - u(i, j))*u.grid.one_over_dy;
					}
					else if (j == j_end)
					{
						grad_x = (u(i + 1, j) - u(i, j))*u.grid.one_over_dx, grad_y = (u(i, j) - u(i, j - 1))*u.grid.one_over_dy;
					}
					else
					{
						grad_x = (u(i + 1, j) - u(i, j))*u.grid.one_over_dx, grad_y = (u(i, j + 1) - u(i, j - 1))*u.grid.one_over_2dy;
					}
				}
				else if (i == i_end)
				{
					if (j == j_start)
					{
						grad_x = (u(i, j) - u(i - 1, j))*u.grid.one_over_dx, grad_y = (u(i, j + 1) - u(i, j))*u.grid.one_over_dy;
					}
					else if (j == j_end)
					{
						grad_x = (u(i, j) - u(i - 1, j))*u.grid.one_over_dx, grad_y = (u(i, j) - u(i, j - 1))*u.grid.one_over_dy;
					}
					else
					{
						grad_x = (u(i, j) - u(i - 1, j))*u.grid.one_over_dx, grad_y = (u(i, j + 1) - u(i, j - 1))*u.grid.one_over_2dy;
					}
				}
				else if (j == j_start)
				{
					grad_x = (u(i + 1, j) - u(i - 1, j))*u.grid.one_over_2dx, grad_y = (u(i, j + 1) - u(i, j))*u.grid.one_over_dy;
				}
				else if (j == j_end)
				{
					grad_x = (u(i + 1, j) - u(i - 1, j))*u.grid.one_over_2dx, grad_y = (u(i, j) - u(i, j - 1))*u.grid.one_over_dy;
				}
				else
				{
					grad_x = (u(i + 1, j) - u(i - 1, j))*u.grid.one_over_2dx, grad_y = (u(i, j + 1) - u(i, j - 1))*u.grid.one_over_2dy;
				}
								
				grad.x = grad_x;
				grad.y = grad_y;
			}
		}
	}

	void ComputeGradient(const FIELD_STRUCTURE_2D<T>& solution, const int& thread_id)
	{
		BEGIN_GRID_ITERATION_2D(solution.partial_grids[thread_id]);
		{
			cout << i << endl;
			cout << j << endl;

			VT& grad(grad_u(i, j));

			T grad_x = (solution(i + 1, j) - solution(i - 1, j))*u.grid.one_over_2dx, grad_y = (solution(i, j + 1) - solution(i, j - 1))*u.grid.one_over_2dy;
			
			grad.x = grad_x;
			grad.y = grad_y;
		}	
		END_GRID_ITERATION_2D;
	}
};