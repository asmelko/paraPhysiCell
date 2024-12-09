#include "custom.h"

#include <cmath>

#include "src/models/morse_position_model.h"
#include "src/original/core/signal_behavior.h"
#include "src/original/modules/geometry.h"
#include "src/random.h"

void create_cell_types(builder& builder)
{
	/*
	   Put any modifications to default cell definition here if you
	   want to have "inherited" by other cell types.

	   This is a good place to set default functions.
	*/

	auto& cell_defaults = builder.get_default_cell_definition();

	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based;


	// This parses the cell definitions in the XML config file and builds the map of cell definitions and summarizes the
	// setup.

	// auto& cell_definitions = builder.get_cell_definitions();

	/*
	   Put any modifications to individual cell definitions here.

	   This is a good place to set custom functions.
	*/
}

void setup_microenvironment(microenvironment_builder&)
{
	// set domain parameters

	// put any custom code to set non-homogeneous initial conditions or
	// extra Dirichlet nodes here.
}

void setup_tissue(environment& e, User_Parameters& parameters, const pugi::xml_node& config_root)
{
	// fill data according to custom parameters

	for (index_t i = 0; i < e.cell_definitions_count; i++)
	{
		index_t N = e.cell_definitions[i]->custom_data["N"];
		real_t cell_radius = e.cell_definitions[i]->custom_data["cell_radius"];

		real_t element_radius;

		if (e.cell_definitions[i]->custom_data["shape"] == 1) // circle
		{
			real_t packing_fraction = 0.9;
			element_radius = cell_radius * std::sqrt(packing_fraction / N);
		}
		else // square
		{
			element_radius = cell_radius / std::sqrt(N);
		}

		e.cell_definitions[i]->phenotype.geometry.radius() = element_radius;

		std::cout << "element_radius: " << element_radius << std::endl;
	}

	for (index_t i = 0; i < e.cell_definitions_count; i++)
	{
		for (index_t j = 0; j < e.cell_definitions_count; j++)
		{
			index_t min = std::min(i, j);
			index_t max = std::max(i, j);

			e.inter_scaling_factors[i * e.cell_definitions_count + j] =
				parameters.doubles("scaling_factor_" + std::to_string(min) + "_" + std::to_string(max));

			e.inter_stiffnesses[i * e.cell_definitions_count + j] =
				parameters.doubles("stiffness_" + std::to_string(min) + "_" + std::to_string(max));

			e.inter_equilibrium_distances[i * e.cell_definitions_count + j] =
				parameters.doubles("equilibrium_distance_multiplier_" + std::to_string(min) + "_" + std::to_string(max))
				* (e.cell_definitions[i]->phenotype.geometry.radius()
				   + e.cell_definitions[j]->phenotype.geometry.radius());
		}
	}

	if (parameters.strings("potential") == "morse")
	{
		e.position = std::make_unique<morse_position_model>();
	}

	e.cell_definitions[1]->functions.custom_cell_rule = [](cell& c) {
		// forced to move in a direction
		c.velocity()[1] -= 1;
	};

	{
		index_t N = e.cell_definitions[0]->custom_data["N"];
		real_t element_radius = e.cell_definitions[0]->phenotype.geometry.radius();


		// e.cell_definitions[0]->phenotype.mechanics.attachment_elastic_constant();

		// create a rectangle of cells
		real_t x_offset_begin = 0;
		real_t y_offset_begin = 100;
		real_t x_offset = x_offset_begin;
		real_t y_offset = y_offset_begin;
		index_t row = 0;
		index_t col = 0;
		for (index_t i = 0; i < N; i++)
		{
			auto cell = e.get_container().create_cell(e.cell_defaults());
			cell->residency() = 1;

			point_t<real_t, 3> position = { 0, 0, 0 };

			position[0] = x_offset;
			position[1] = y_offset;
			position[2] = 0;

			cell->assign_position(position);

			std::cout << "placing cell " << i << " at position " << position[0] << " " << position[1] << std::endl;
			x_offset += std::sqrt(3) * element_radius;
			y_offset += std::pow(-1, col) * element_radius;
			col++;
			if ((i + 1) % (int)std::sqrt(N) == 0)
			{
				row++;
				col = 0;
				x_offset = x_offset_begin;
				y_offset = y_offset_begin + 2. * row * element_radius;
			}
		}
	}

	{
		index_t N = e.cell_definitions[0]->custom_data["N"];
		real_t element_radius = e.cell_definitions[0]->phenotype.geometry.radius();


		// e.cell_definitions[0]->phenotype.mechanics.attachment_elastic_constant();

		// create a rectangle of cells
		real_t x_offset_begin = 170;
		real_t y_offset_begin = 100;
		real_t x_offset = x_offset_begin;
		real_t y_offset = y_offset_begin;
		index_t row = 0;
		index_t col = 0;
		for (index_t i = 0; i < N; i++)
		{
			auto cell = e.get_container().create_cell(e.cell_defaults());
			cell->residency() = 2;

			point_t<real_t, 3> position = { 0, 0, 0 };

			position[0] = x_offset;
			position[1] = y_offset;
			position[2] = 0;

			cell->assign_position(position);

			std::cout << "placing cell " << i << " at position " << position[0] << " " << position[1] << std::endl;
			x_offset += std::sqrt(3) * element_radius;
			y_offset += std::pow(-1, col) * element_radius;
			col++;
			if ((i + 1) % (int)std::sqrt(N) == 0)
			{
				row++;
				col = 0;
				x_offset = x_offset_begin;
				y_offset = y_offset_begin + 2. * row * element_radius;
			}
		}
	}

	{
		index_t N = e.cell_definitions[1]->custom_data["N"];
		real_t element_radius = e.cell_definitions[1]->phenotype.geometry.radius();


		// create a rectangle of cells
		real_t x_offset = 160 - element_radius;
		real_t y_offset = 400;
		for (index_t i = 0; i < N; i++)
		{
			auto cell = e.get_container().create_cell(*e.cell_definitions[1]);
			cell->residency() = 3;

			point_t<real_t, 3> position = { 0, 0, 0 };

			position[0] = x_offset;
			position[1] = y_offset;
			position[2] = 0;

			cell->assign_position(position);

			std::cout << "placing cell " << i << " at position " << position[0] << " " << position[1] << std::endl;
			x_offset += 2. * element_radius;
			if ((i + 1) % (int)(N / std::sqrt(N)) == 0 && i != 0)
			{
				x_offset = 160 - element_radius;
				y_offset += 2. * element_radius;
			}
		}
	}
}

cell_coloring_funct_t get_my_coloring_function(User_Parameters&) { return paint_by_number_cell_coloring; }
