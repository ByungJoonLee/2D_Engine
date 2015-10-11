#pragma once

#include "OPENGL_OBJECT_BASE.h"

class OPENGL_LEVELSET : public OPENGL_OBJECT_BASE
{
public: // Essential Data
	bool								is_polygonized;
	int									num_threads;
	float								grid_scale;
	MULTITHREADING*						multithreading;
	MARCHING_SQUARES_ALGORITHM			levelset_polygonizer, levelset_polygonizer_2;
	
	int									length_of_drawing_object;
	T									min_level, max_level;

	ARRAY<MARCHING_SQUARES_ALGORITHM>	levelset_polygonizer_array;
	LEVELSET_2D*						levelset_object;
	LEVELSET_2D*						second_levelset_object;

public: // Enumerate
	enum SLEVELSET_DRAW_TYPE
	{
		LEVELSET_DRAW_HIDE			= 0x0000,
		LEVELSET_DRAW_SOLID			= 0x0001,
		LEVELSET_DRAW_WIRE			= 0x0002,
		LEVELSET_DRAW_WIRE_SOLID	= LEVELSET_DRAW_SOLID | LEVELSET_DRAW_WIRE,
		LEVELSET_DRAW_WIRE_DOUBLE	= 0x0004,
		LEVELSET_DRAW_WIRE_ARRAY    = 0x0005,
	};

public: // Constructor and Destructor
	OPENGL_LEVELSET(const char* display_name, OPENGL_DRIVER* driver, LEVELSET_2D* levelset_object_input, MULTITHREADING* multithreading_input, float grid_scale_input)
		: OPENGL_OBJECT_BASE(display_name, driver), levelset_object(levelset_object_input), second_levelset_object(0), num_threads(multithreading_input->num_threads), multithreading(multithreading_input), grid_scale(grid_scale_input), is_polygonized(false), length_of_drawing_object((int)1),
		  min_level((T)0), max_level((T)0)
	{
		RegisterDrawType((int) OPENGL_LEVELSET::LEVELSET_DRAW_HIDE, "HIDE");
		RegisterDrawType((int) OPENGL_LEVELSET::LEVELSET_DRAW_SOLID, "SOLID");
		RegisterDrawType((int) OPENGL_LEVELSET::LEVELSET_DRAW_WIRE, "WIRE");
		RegisterDrawType((int) OPENGL_LEVELSET::LEVELSET_DRAW_WIRE_SOLID, "WIRE_SOLID");
		RegisterDrawType((int) OPENGL_LEVELSET::LEVELSET_DRAW_WIRE_DOUBLE, "WIRE_DOUBLE");
		RegisterDrawType((int) OPENGL_LEVELSET::LEVELSET_DRAW_WIRE_ARRAY, "WIRE_ARRAY");
		SetDrawType((int) LEVELSET_DRAW_WIRE);
		is_levelset = true;
	}

	OPENGL_LEVELSET(const char* display_name, OPENGL_DRIVER* driver, LEVELSET_2D* levelset_object_input, LEVELSET_2D* second_levelset_object_input, MULTITHREADING* multithreading_input, float grid_scale_input)
		: OPENGL_OBJECT_BASE(display_name, driver), levelset_object(levelset_object_input), second_levelset_object(second_levelset_object_input), num_threads(multithreading_input->num_threads), multithreading(multithreading_input), grid_scale(grid_scale_input), is_polygonized(false), length_of_drawing_object((int)1),
		  min_level((T)0), max_level((T)0)
	{
		RegisterDrawType((int) OPENGL_LEVELSET::LEVELSET_DRAW_HIDE, "HIDE");
		RegisterDrawType((int) OPENGL_LEVELSET::LEVELSET_DRAW_SOLID, "SOLID");
		RegisterDrawType((int) OPENGL_LEVELSET::LEVELSET_DRAW_WIRE, "WIRE");
		RegisterDrawType((int) OPENGL_LEVELSET::LEVELSET_DRAW_WIRE_SOLID, "WIRE_SOLID");
		RegisterDrawType((int) OPENGL_LEVELSET::LEVELSET_DRAW_WIRE_DOUBLE, "WIRE_DOUBLE");
		RegisterDrawType((int) OPENGL_LEVELSET::LEVELSET_DRAW_WIRE_ARRAY, "WIRE_ARRAY");
		SetDrawType((int) LEVELSET_DRAW_WIRE);
		is_levelset = true;
	}

	~OPENGL_LEVELSET(void)
	{}

public: // Member Functions
	virtual int NextDrawType()
	{
		switch(GetDrawType())
		{
		case LEVELSET_DRAW_HIDE:
			SetDrawType((int) LEVELSET_DRAW_SOLID);
			break;
		case LEVELSET_DRAW_SOLID:
			SetDrawType((int) LEVELSET_DRAW_WIRE);
			break;
		case LEVELSET_DRAW_WIRE:
			SetDrawType((int) LEVELSET_DRAW_WIRE_SOLID);
			break;
		case LEVELSET_DRAW_WIRE_DOUBLE:
			SetDrawType((int) LEVELSET_DRAW_WIRE_DOUBLE);
			break;
		case LEVELSET_DRAW_WIRE_ARRAY:
			SetDrawType((int) LEVELSET_DRAW_WIRE_ARRAY);
			break;
		case LEVELSET_DRAW_WIRE_SOLID:
			SetDrawType((int) LEVELSET_DRAW_HIDE);
			break;
		}

		return GetDrawType();
	}

	virtual int PreviousDrawType()
	{
		switch(GetDrawType())
		{
		case LEVELSET_DRAW_HIDE:
			SetDrawType((int) LEVELSET_DRAW_WIRE_SOLID);
			break;
		case LEVELSET_DRAW_SOLID:
			SetDrawType((int) LEVELSET_DRAW_HIDE);
			break;
		case LEVELSET_DRAW_WIRE:
			SetDrawType((int) LEVELSET_DRAW_SOLID);
			break;
		case LEVELSET_DRAW_WIRE_DOUBLE:
			SetDrawType((int) LEVELSET_DRAW_WIRE_DOUBLE);
			break;
		case LEVELSET_DRAW_WIRE_ARRAY:
			SetDrawType((int) LEVELSET_DRAW_WIRE_ARRAY);
			break;
		case LEVELSET_DRAW_WIRE_SOLID:
			SetDrawType((int) LEVELSET_DRAW_WIRE);
			break;
		}

		return GetDrawType();
	}

	virtual void Update()
	{
		if (!levelset_object)
		{
			return;
		}
				
		if (GetDrawType() == LEVELSET_DRAW_HIDE)
		{
			is_polygonized = false;
		}
		else if (GetDrawType() == LEVELSET_DRAW_WIRE_DOUBLE)
		{
			levelset_polygonizer.Initialize(multithreading, levelset_object->grid, grid_scale);
			levelset_polygonizer_2.Initialize(multithreading, levelset_object->grid, grid_scale);

			LEVELSET_2D drawing_object, drawing_object_2;
			drawing_object.Initialize(levelset_object->grid, 2, multithreading);
			drawing_object_2.Initialize(levelset_object->grid, 2, multithreading);
			
			// Drawing level choosing option
			GRID_ITERATION_2D(drawing_object.grid)
			{
				drawing_object(i, j) = 	levelset_object->arr(i, j) - levelset_object->epsilon;
				drawing_object_2(i, j) = levelset_object->arr(i, j) + levelset_object->epsilon;
			}
			
			levelset_polygonizer.Polygonize(drawing_object, true);
			levelset_polygonizer_2.Polygonize(drawing_object_2, true);

			//levelset_polygonizer.Polygonize(*(levelset_object), true);

			is_polygonized = true;
		}
		else if (GetDrawType() == LEVELSET_DRAW_WIRE_ARRAY)
		{
			levelset_polygonizer_array.Initialize(2*length_of_drawing_object);
			for (int i = 0; i < 2*length_of_drawing_object; i++)
			{
				levelset_polygonizer_array[i].Initialize(multithreading, levelset_object->grid, grid_scale);
			}
			
			/*levelset_polygonizer.Initialize(multithreading, levelset_object->grid, grid_scale);
			levelset_polygonizer_2.Initialize(multithreading, levelset_object->grid, grid_scale);*/
						
			ARRAY<LEVELSET_2D> drawing_object_array;
			drawing_object_array.Initialize(2*length_of_drawing_object);
			for (int i = 0; i < drawing_object_array.length; i++)
			{
				drawing_object_array[i].Initialize(levelset_object->grid, 2, multithreading);
			}

			for (int k = 0; k < 2*length_of_drawing_object; k++)
			{
				if (k < length_of_drawing_object)
				{
					GRID_ITERATION_2D(drawing_object_array[k].grid)
					{
						(drawing_object_array[k])(i, j) = levelset_object->arr(i, j) - k*(max_level - min_level)/length_of_drawing_object;
					}
				}
				else
				{
					GRID_ITERATION_2D(drawing_object_array[k].grid)
					{
						(drawing_object_array[k])(i, j) = second_levelset_object->arr(i, j) + (k - length_of_drawing_object)*(max_level - min_level)/length_of_drawing_object;
					}
				}
				
			}

			for (int i = 0; i < 2*length_of_drawing_object; i++)
			{
				levelset_polygonizer_array[i].Polygonize(drawing_object_array[i], true);
			}
			
			is_polygonized = true;
		}
		else
		{
			levelset_polygonizer.Initialize(multithreading, levelset_object->grid, grid_scale);
			
			LEVELSET_2D drawing_object;
			drawing_object.Initialize(levelset_object->grid, 2, multithreading);
			
			// Drawing level choosing option
			GRID_ITERATION_2D(drawing_object.grid)
			{
				drawing_object(i, j) = 	levelset_object->arr(i, j) - levelset_object->epsilon;
			}
			
			levelset_polygonizer.Polygonize(drawing_object, true);

			//levelset_polygonizer.Polygonize(*(levelset_object), true);

			is_polygonized = true;
		}
	}

	virtual void Render()
	{
		if (!levelset_object)
		{
			return;
		}

		if (GetDrawType() == LEVELSET_DRAW_SOLID)
		{
			RenderSolid();
		}

		if (GetDrawType() == LEVELSET_DRAW_WIRE)
		{
			RenderWire();
		}

		if (GetDrawType() == LEVELSET_DRAW_WIRE_DOUBLE)
		{
			RenderWireDouble();
		}

		if (GetDrawType() == LEVELSET_DRAW_WIRE_ARRAY)
		{
			RenderWireArray();
		}
	}

	void RenderSolid()
	{
		if (!is_polygonized)
		{
			Update();
		}

		GetDriver()->SetRenderStatesByMaterial(material);
		for (int thread_id = 0; thread_id < num_threads; thread_id++)
		{
			if (levelset_polygonizer.surfaces[thread_id])
			{
				levelset_polygonizer.surfaces[thread_id]->DrawTriangles(true);
			}
		}
	}

	void RenderWire()
	{
		if (!is_polygonized)
		{
			return;
		}

		GetDriver()->SetDefaultRenderStatesLineMode();
		for (int thread_id = 0; thread_id < num_threads; thread_id++)
		{
			if (levelset_polygonizer.surfaces[thread_id])
			{
				levelset_polygonizer.surfaces[thread_id]->DrawEdges();
				//levelset_polygonizer.surfaces[thread_id]->DrawVertices();
			}
		}
	}

	void RenderWireDouble()
	{
		if (!is_polygonized)
		{
			return;
		}

		GetDriver()->SetDefaultRenderStatesLineMode();
		for (int thread_id = 0; thread_id < num_threads; thread_id++)
		{
			if (levelset_polygonizer.surfaces[thread_id])
			{
				levelset_polygonizer.surfaces[thread_id]->DrawEdges();
				//levelset_polygonizer.surfaces[thread_id]->DrawVertices();
			}
			
			if (levelset_polygonizer_2.surfaces[thread_id])
			{
				levelset_polygonizer_2.surfaces[thread_id]->DrawEdges();
			}
		}
	}

	void RenderWireArray()
	{
		if (!is_polygonized)
		{
			return;
		}

		GetDriver()->SetDefaultRenderStatesLineMode();
		for (int thread_id = 0; thread_id < num_threads; thread_id++)
		{
			for (int i = 0; i < 2*length_of_drawing_object; i++)
			{
				if (levelset_polygonizer_array[i].surfaces[thread_id])
				{
					levelset_polygonizer_array[i].surfaces[thread_id]->DrawEdges();
				}
			}
		}
	}

	void SetLengthForDrawingObject(const int& length_of_drawing_object_input)
	{
		length_of_drawing_object = length_of_drawing_object_input;
	}

	void SetMinMaxLevelForDrawing(const T& min_level_input, const T& max_level_input)
	{
		min_level = min_level_input;
		max_level = max_level_input;
	}
	
	MARCHING_SQUARES_ALGORITHM* MarchingSquareAlgorithm()
	{
		return& levelset_polygonizer;
	}
};