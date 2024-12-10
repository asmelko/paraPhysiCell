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

void make_packed_square(environment& e, cell_definition* cd, index_t residency, point_t<real_t, 3> starting_position)
{
	index_t N = cd->custom_data["N"];
	real_t element_radius = cd->phenotype.geometry.radius();

	auto position = starting_position;

	index_t row = 0;
	index_t col = 0;
	for (index_t i = 0; i < N; i++)
	{
		auto cell = e.get_container().create_cell(*cd);
		cell->residency() = residency;

		cell->assign_position(position);

		std::cout << "placing cell " << i << " at position " << position[0] << " " << position[1] << std::endl;
		position[0] += std::sqrt(3) * element_radius;
		position[1] += std::pow(-1, col) * element_radius;
		col++;
		if ((i + 1) % (int)std::sqrt(N) == 0)
		{
			row++;
			col = 0;
			position[0] = starting_position[0];
			position[1] = starting_position[1] + 2. * row * element_radius;
		}
	}
}

void make_circle(environment& e, cell_definition* cd, index_t residency, point_t<real_t, 3> starting_position)
{
	index_t N = cd->custom_data["N"];
	real_t element_radius = cd->phenotype.geometry.radius();
	real_t cell_radius = cd->custom_data["cell_radius"];

	auto position = starting_position;

	index_t row = 0;

	while (true)
	{
		index_t line = std::sqrt(cell_radius * cell_radius - row * row * 4 * element_radius * element_radius);
		index_t line_count = std::round(line / element_radius);
		index_t offset = cell_radius - (line_count * element_radius);

		if (offset < 0)
		{
			break;
		}

		position[0] = starting_position[0] + offset;

		for (index_t i = 0; i < line_count; i++)
		{
			for (auto pos : (row == 0 ? std::vector<index_t> { 1 } : std::vector<index_t> { -1, 1 }))
			{
				auto cell = e.get_container().create_cell(*cd);
				cell->residency() = residency;

				position[1] = starting_position[1] + pos * row * 2 * element_radius;

				cell->assign_position(position);
				std::cout << "placing cell " << i << " at position " << position[0] << " " << position[1] << std::endl;

				if (--N == 0)
				{
					return;
				}
			}

			position[0] += 2. * element_radius;
		}

		row++;
	}
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
			real_t packing_fraction = 0.5;
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

	make_circle(e, e.cell_definitions[0].get(), 0, { -50, 0 });
	make_circle(e, e.cell_definitions[0].get(), 1, { 100, 0 });
	make_packed_square(e, e.cell_definitions[0].get(), 2, { 300, 0 });
	make_packed_square(e, e.cell_definitions[0].get(), 3, { 450, 0 });

	make_circle(e, e.cell_definitions[1].get(), 4, { -60, 200 });
	make_circle(e, e.cell_definitions[1].get(), 5, { 110, 200 });
	make_packed_square(e, e.cell_definitions[1].get(), 6, { 280, 200 });
	make_packed_square(e, e.cell_definitions[1].get(), 7, { 460, 200 });
}

cell_coloring_funct_t get_my_coloring_function(User_Parameters&) { return paint_by_number_cell_coloring; }
