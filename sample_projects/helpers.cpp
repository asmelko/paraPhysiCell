#include "helpers.h"

#include <cmath>

#include "src/models/kelvin_voigt_model.h"
#include "src/models/morse_kv_membrane_model.h"
#include "src/models/morse_position_model.h"

using namespace biofvm;
using namespace physicell;

void physicell::make_packed_square(environment& e, cell_definition* cd, index_t residency,
								   point_t<real_t, 3> starting_position)
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

void physicell::make_circle(environment& e, cell_definition* cd, index_t residency,
							point_t<real_t, 3> starting_position)
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

				position[1] = starting_position[1] + pos * row * 1.7 * element_radius;

				cell->assign_position(position);
				std::cout << "placing cell " << i << " at position " << position[0] << " " << position[1] << std::endl;

				if (--N == 0)
				{
					return;
				}
			}

			position[0] += 1.99 * element_radius;
		}

		row++;
	}
}

void physicell::setup_potential_parameters(environment& e, User_Parameters& parameters)
{
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

		e.membrane_stretching_factors[i] = e.cell_definitions[i]->custom_data["membrane_stretching_factor"];
		e.membrane_stretched_spring_stepup[i] =
			e.cell_definitions[i]->custom_data["membrane_stretched_spring_stepup"];
	}

	if (parameters.strings("potential") == "morse")
	{
		e.position = std::make_unique<morse_position_model>();
	}
	else if (parameters.strings("potential") == "kelvin_voigt")
	{
		e.position = std::make_unique<kelvin_voigt_model>();
	}
	else if (parameters.strings("potential") == "morse_kv_membrane")
	{
		e.position = std::make_unique<morse_kv_membrane_position_model>();
	}
}
