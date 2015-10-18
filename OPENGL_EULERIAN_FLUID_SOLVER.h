#pragma once

#include "OPENGL_SOLVER_BASE.h"
#include "EULERIAN_FLUID_SOLVER_2D.h"
#include "OPENGL_LEVELSET.h"
#include "OPENGL_SCALARFIELD.h"
#include "OPENGL_VECTORFIELD.h"
#include <vector>

class EULERIAN_FLUID_SOLVER_2D;
class OPENGL_LEVELSET;
class OPENGL_SCALARFIELD;

class OPENGL_EULERIAN_FLUID_SOLVER : public OPENGL_SOLVER_BASE
{
public: // Essential Data
	OPENGL_DRIVER*					driver;
	EULERIAN_FLUID_SOLVER_2D*		eulerian_solver;

	OPENGL_LEVELSET*				water_levelset;
	//OPENGL_LEVELSET*				vortex_levelset;

	OPENGL_SCALARFIELD*				vortex_levelset;
	OPENGL_SCALARFIELD*				pressure_field;
	OPENGL_VECTORFIELD*				velocity_field;
	

public: // Constructor and Destructor
	OPENGL_EULERIAN_FLUID_SOLVER(OPENGL_DRIVER* driver_input, EULERIAN_FLUID_SOLVER_2D* eulerian_object, MULTITHREADING* multithreading_input, float grid_scale)
		: driver(driver_input), eulerian_solver(eulerian_object), water_levelset(0), vortex_levelset(0), velocity_field(0), pressure_field(0)
	{
		if (eulerian_solver->vortex_sheet_problem == true && eulerian_solver->test_number == 3)
		{
			water_levelset = new OPENGL_LEVELSET("EULERIAN_WATER_LEVELSET", driver, eulerian_solver->water_levelset, eulerian_solver->second_levelset, multithreading_input, grid_scale);
			AddObject(water_levelset);
		}
		else if (eulerian_solver->bifurcation_test == true)
		{
			vortex_levelset = new OPENGL_SCALARFIELD("VORTEX_FUNCTION", driver, &eulerian_solver->vortex_levelset->signed_distance_field);
			AddObject(vortex_levelset);
		}
		else
		{
			water_levelset = new OPENGL_LEVELSET("EULERIAN_WATER_LEVELSET", driver, eulerian_solver->water_levelset, multithreading_input, grid_scale);
			AddObject(water_levelset);
		}
				
		if (eulerian_solver->vortex_sheet_problem)
		{
			if (eulerian_solver->test_number == 2)
			{
				water_levelset->SetDrawType(OPENGL_LEVELSET::LEVELSET_DRAW_WIRE_DOUBLE);
			}
			if (eulerian_solver->test_number == 3)
			{
				water_levelset->SetDrawType(OPENGL_LEVELSET::LEVELSET_DRAW_WIRE_ARRAY);
				water_levelset->SetLengthForDrawingObject(10);
				water_levelset->SetMinMaxLevelForDrawing(0, 0.04);
			}
		}

		/*if (eulerian_solver->vortex_sheet_problem)
		{
			vortex_levelset = new OPENGL_LEVELSET("VORTEX_WATER_LEVELSET", driver, eulerian_solver->vortex_levelset, multithreading_input, grid_scale);
			AddObject(vortex_levelset);
		}
		else
		{
			water_levelset = new OPENGL_LEVELSET("EULERIAN_WATER_LEVELSET", driver, eulerian_solver->water_levelset, multithreading_input, grid_scale);
			AddObject(water_levelset);
		}*/
		
		velocity_field = new OPENGL_VECTORFIELD("VELOCITY_FIELD", driver, eulerian_solver->water_velocity_field, eulerian_solver->water_velocity_field_mac_x, eulerian_solver->water_velocity_field_mac_y);
		AddObject(velocity_field);
		if (eulerian_solver->vortex_sheet_problem || eulerian_solver->bifurcation_test)
		{
			velocity_field->SetDrawType(OPENGL_VECTORFIELD::VECTORFIELD_DRAW_SHOW);
		}
		else
		{
			velocity_field->SetDrawType(OPENGL_VECTORFIELD::VECTORFIELD_DRAW_SHOW_MAC);
		}
		
		if (eulerian_solver->bifurcation_test)
		{
			pressure_field = new OPENGL_SCALARFIELD("STREAM_FUNCTION", driver, &eulerian_solver->water_projection->pressure_field);
			AddObject(pressure_field);
			pressure_field->SetDrawType(OPENGL_SCALARFIELD::SCALARFIELD_DRAW_SHOW);
		}
		else
		{
			pressure_field = new OPENGL_SCALARFIELD("PRESSURE_FIELD", driver, &eulerian_solver->water_projection->pressure_field);
			AddObject(pressure_field);
			pressure_field->SetDrawType(OPENGL_SCALARFIELD::SCALARFIELD_DRAW_SHOW);
		}
	}

	~OPENGL_EULERIAN_FLUID_SOLVER(void)
	{}

public: // Member Functions
	virtual void LoadOpenGLSettings(SCRIPT_BLOCK* script_block)
	{
		// Water Levelset
		SCRIPT_BLOCK* solver_block = script_block->SearchBlock("OPENGL_EULERIAN_FLUID_SOLVER");
		if (!solver_block)
		{
			return;
		}

		SCRIPT_BLOCK* water_levelset_block = solver_block->SearchBlock("WATER_LEVELSET");
		if (water_levelset && water_levelset_block)
		{
			OPENGL_MATERIAL material;
			const char* solid_color_name = water_levelset_block->GetString("solid_material_name", 0);
			if (solid_color_name)
			{
				material.SetMaterialByName(solid_color_name);
			}
			else
			{
				VECTOR_3D<T> ambient = water_levelset_block->GetVector3("solid_ambient_color");
				VECTOR_3D<T> diffuse = water_levelset_block->GetVector3("solid_diffuse_color");
				VECTOR_3D<T> specular = water_levelset_block->GetVector3("solid_specular_color");
				GLfloat shininess = water_levelset_block->GetFloat("solid_shininess");

				material.ambient_color.Set(ambient.x, ambient.y, ambient.z);
				material.diffuse_color.Set(diffuse.x, diffuse.y, diffuse.z);
				material.specular_color.Set(specular.x, specular.y, specular.z);
				material.shininess = shininess;
			}

			water_levelset->SetMaterial(material);

			// Display given variables
			cout << "-------------- OPENGL EULERIAN FLUID SOLVER - LEVELSET --------------" << endl;
			cout << "solid material name : " << solid_color_name << endl;
		}
	}
};