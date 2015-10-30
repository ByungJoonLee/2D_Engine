#pragma once

#include "COMMON_DEFINITION.h"
#include "FIELD_STRUCTURE_2D.h"
#include "LEVELSET_2D.h"
#include "GS_COEFFICIENTS.h"
#include "POISSON_SOLVER.h"

class MULTILEVEL_SOLVER
{
public: 
	MULTITHREADING			*multithreading;
	
	int						projection_type;
	
public: // Constructor and Destructor
	MULTILEVEL_SOLVER(MULTITHREADING* multithreading_input);
	~MULTILEVEL_SOLVER()
	{}
};