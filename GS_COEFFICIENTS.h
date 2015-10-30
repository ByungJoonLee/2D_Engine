#pragma once

#include "COMMON_DEFINITION.h"

class GS_COEFFICIENTS
{
public: // Constructor and Destructor
	GS_COEFFICIENTS()
	{
		Initialize();
	}

public: 
	T plus[2], minus[2];
	T center, inv_center;

public: 
	inline void Initialize()
	{
		plus[0] = (T)0;
		plus[1] = (T)0;
		minus[0] = (T)0;
		minus[1] = (T)0;
		center = (T)0;
		inv_center = (T)1;
	}

	inline void operator = (const GS_COEFFICIENTS& v)
	{
		plus[0] = v.plus[0];
		plus[1] = v.plus[1];
		minus[0] = v.minus[0];
		minus[1] = v.minus[1];
		center = (T)0;
		inv_center = (T)1;
	}

	void Cout()
	{
		std::cout << center << " " << inv_center << " " << plus[0] << " " << plus[1] << " " << " " << minus[0] << " " << minus[1] << std::endl;
	}
};
