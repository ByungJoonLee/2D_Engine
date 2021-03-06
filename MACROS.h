#pragma once

enum  POISSON_SOLVER_TYPE												{NO_SOLVER, CG, PCG, GS, BICG, SOR, MULTIGRID};

#define PI																(T)3.141592
#define BC_FULL															0
#define BC_DIR															-1
#define BC_OBJ															-2
#define BC_NULL															-3
#define BC_NEUM															-4
#define BC_IMPLICIT														-5
#define BC_PER															-6

#define CLAMP(v, min, max)												((v) < (min) ? (min) : ((v) > (max) ? (max) : (v)))

#define POW2(a)															((a)*(a))
#define POW3(a)															((a)*(a)*(a))

#define MAX(a, b)														((a) > (b) ? (a) : (b))
#define MIN(a, b)														((a) > (b) ? (b) : (a))

#define INCREASING_SORT2(a, b, a1, a2)									if (a <= b){a1 = a; a2 = b;}                 \
																		else{a1 = b; a2 = a;}                        

#define INCREASING_SORT3(a, b, c, a1, a2, a3)							if(a <= b){											\
																			if(b <= c) {a1 = a; a2 = b; a3 = c;}			\
																			else if(a <= c) {a1 = a; a2 = c; a3 = b;}		\
																			else {a1 = c; a2 = a; a3 = b;}}					\
																		else{												\
																			if(a <= c) {a1 = b; a2 = a; a3 = c;}			\
																			else if(b <= c) {a1 = b; a2 = c; a3 = a;}		\
																			else {a1 = c; a3 = b; a3 = a;}}

#define DELETE_POINTER(pointer)											if (pointer != 0) {delete pointer; pointer = 0;}
#define DELETE_ARRAY(pointer)											if (pointer != 0) {delete [] pointer; pointer = 0;}

#define LOOPS_1D(i, i_start, i_end)										for((i) = (i_start); (i) <= (i_end); ++(i))

#define LOOPS_2D(i, j, i_start, j_start, i_end, j_end)					for((j) = (j_start) ; (j) <= (j_end); ++(j)) for((i) = (i_start); (i) <= (i_end); ++(i))

#define GRID_ITERATION_1D(grid_1d_input)								int i_start = (grid_1d_input).i_start, i_end = (grid_1d_input).i_end, i; for(i = i_start; i <= i_end; ++i)

#define GRID_ITERATION_2D(grid_2d_input)								for(int j_start = (grid_2d_input).j_start, j_end = (grid_2d_input).j_end, i_start = (grid_2d_input).i_start, i_end = (grid_2d_input).i_end, i, j = j_start; j <= j_end; ++j) for(i = i_start; i <= i_end; ++i)

#define BEGIN_GRID_ITERATION_1D(grid_1d_input)							{GRID_STRUCTURE_1D& grid_1d_itr(grid_1d_input);																								\
																		int i(0);																																	\
																		const int i_start = grid_1d_itr.i_start, i_end = grid_1d_itr.i_end;				\
																		for (int i = i_start; i <= i_end; ++i) 

#define BEGIN_GRID_ITERATION_2D(grid_2d_input)							{GRID_STRUCTURE_2D& grid_2d_itr(grid_2d_input);																								\
																		int i(0), j(0);																																	\
																		const int j_start = grid_2d_itr.j_start, j_end = grid_2d_itr.j_end, i_start = grid_2d_itr.i_start, i_end = grid_2d_itr.i_end;				\
																		for (int j = j_start; j <= j_end; ++j) for (int i = i_start; i <= i_end; ++i) 

#define BEGIN_GRID_ITERATION_2D_ARR(grid_2d_input, arr)					{GRID_STRUCTURE_2D &grid_2d_itr(grid_2d_input);	\
																		int i(0), j(0), arr_ix(0);						\
																		const int j_start = grid_2d_itr.j_start, j_end = grid_2d_itr.j_end, i_start = grid_2d_itr.i_start, i_end = grid_2d_itr.i_end;	\
																		for(j = j_start; j <= j_end; ++j) for(i = i_start, arr_ix=arr.Index1D(i, j); i <= i_end; ++i, ++arr_ix)

#define BEGIN_GRID_ITERATION_2D_ARR2(grid_2d_input, arr1, arr2)			{GRID_STRUCTURE_2D &grid_itr(grid_2d_input);	\
																		int i, j, arr1_ix, arr2_ix;						\
																		const int j_start = grid_itr.j_start, j_end = grid_itr.j_end, i_start = grid_itr.i_start, i_end = grid_itr.i_end;	\
																		for(j = j_start; j <= j_end; ++j) for(i = i_start, arr1_ix = arr1.Index1D(i, j), arr2_ix = arr2.Index1D(i, j); i <= i_end; ++i, ++arr1_ix, ++arr2_ix)

#define BEGIN_GRID_ITERATION_2D_ARR3(grid_2d_input, arr1, arr2, arr3)	{GRID_STRUCTURE_2D &grid_itr(grid_2d_input);	\
																		int i, j, arr1_ix, arr2_ix, arr3_ix;						\
																		const int j_start = grid_itr.j_start, j_end = grid_itr.j_end, i_start = grid_itr.i_start, i_end = grid_itr.i_end;	\
																		for(j = j_start; j <= j_end; ++j) for(i = i_start, arr1_ix = arr1.Index1D(i, j), arr2_ix = arr2.Index1D(i, j), arr3_ix = arr3.Index1D(i ,j); i <= i_end; ++i, ++arr1_ix, ++arr2_ix, ++arr3_ix)

#define END_GRID_ITERATION_1D											multithreading->Sync(thread_id);}

#define END_GRID_ITERATION_2D											multithreading->Sync(thread_id);}

#define BEGIN_HEAD_THREAD_WORK											if(thread_id == 0)
#define END_HEAD_THREAD_WORK											multithreading->Sync(thread_id);

#define HEAD_THREAD_WORK(expression)									if(thread_id == 0){expression;}; multithreading->Sync(thread_id);

#define END_GRID_ITERATION_SUM(sync_value)								multithreading->SyncSum(thread_id, sync_value);}

#define END_GRID_ITERATION_MAX_2D(sync_value)							multithreading->SyncMax(thread_id, sync_value);}

#define PREPARE_FOR_1D_ITERATION(num)									multithreading->SplitDomainIndex1D(thread_id, 0, num);

#define BEGIN_1D_ITERATION												{const int p_start(multithreading->start_ix_1D[thread_id]), p_end(multithreading->end_ix_1D[thread_id]);						\
																		for(int p = p_start; p <= p_end; p++)

#define END_1D_ITERATION												multithreading->Sync(thread_id);}