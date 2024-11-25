#include "custom.h"

#include <cmath>

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
	{
		index_t N = e.cell_defaults().phenotype.mechanics.detachment_rate();
		real_t cell_radius = 100;

		constexpr real_t packing_fraction = 0.7;
		real_t element_radius = cell_radius * std::sqrt(packing_fraction / N);


		e.cell_defaults().phenotype.mechanics.cell_cell_repulsion_strength() = 2.0 * element_radius;
		e.cell_defaults().phenotype.geometry.radius() = element_radius;


		std::cout << "element_radius: " << element_radius << std::endl;
	}

	{
		index_t N = e.cell_definitions[1]->phenotype.mechanics.detachment_rate();
		real_t cell_radius = 50;

		constexpr real_t packing_fraction = 0.7;
		real_t element_radius = cell_radius * std::sqrt(packing_fraction / N);


		e.cell_definitions[1]->phenotype.mechanics.cell_cell_repulsion_strength() = 2.0 * element_radius;
		e.cell_definitions[1]->phenotype.geometry.radius() = element_radius;


		std::cout << "element_radius: " << element_radius << std::endl;
	}


	for (index_t i = 0; i < e.cell_definitions_count; i++)
	{
		for (index_t j = 0; j < e.cell_definitions_count; j++)
		{
			e.cell_definitions[i]->phenotype.mechanics.cell_adhesion_affinities()[j]

				*= (e.cell_definitions[i]->phenotype.geometry.radius()
					+ e.cell_definitions[j]->phenotype.geometry.radius());
		}
	}

	// {
	// 	index_t N = e.cell_defaults().phenotype.mechanics.detachment_rate();
	// 	real_t element_radius = e.cell_defaults().phenotype.geometry.radius();

	// 	// e.cell_defaults().phenotype.mechanics.attachment_elastic_constant();

	// 	// create a rectangle of cells
	// 	real_t x_offset = -200;
	// 	real_t y_offset = -200;
	// 	for (index_t i = 0; i < N; i++)
	// 	{
	// 		auto cell = e.get_container().create_cell(e.cell_defaults());
	// 		cell->residency() = 0;

	// 		point_t<real_t, 3> position = { 0, 0, 0 };

	// 		position[0] = x_offset;
	// 		position[1] = y_offset;
	// 		position[2] = 0;

	// 		cell->assign_position(position);

	// 		std::cout << "placing cell " << i << " at position " << position[0] << " " << position[1] << std::endl;
	// 		x_offset += 2. * element_radius;
	// 		if ((i + 1) % (int)(N / std::sqrt(N)) == 0 && i != 0)
	// 		{
	// 			x_offset = -200;
	// 			y_offset += 2. * element_radius;
	// 		}
	// 	}
	// }

	{
		index_t N = e.cell_defaults().phenotype.mechanics.detachment_rate();
		real_t element_radius = e.cell_defaults().phenotype.geometry.radius();


		// e.cell_defaults().phenotype.mechanics.attachment_elastic_constant();

		// create a rectangle of cells
		real_t x_offset = -20;
		real_t y_offset = -200;
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
			x_offset += 2. * element_radius;
			if ((i + 1) % (int)(N / std::sqrt(N)) == 0 && i != 0)
			{
				x_offset = -20;
				y_offset += 2. * element_radius;
			}
		}
	}

	{
		index_t N = e.cell_defaults().phenotype.mechanics.detachment_rate();
		real_t element_radius = e.cell_defaults().phenotype.geometry.radius();


		// e.cell_defaults().phenotype.mechanics.attachment_elastic_constant();

		// create a rectangle of cells
		real_t x_offset = 150;
		real_t y_offset = -200;
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
			x_offset += 2. * element_radius;
			if ((i + 1) % (int)(N / std::sqrt(N)) == 0 && i != 0)
			{
				x_offset = 150;
				y_offset += 2. * element_radius;
			}
		}
	}

	{
		index_t N = e.cell_definitions[1]->phenotype.mechanics.detachment_rate();
		real_t element_radius = e.cell_definitions[1]->phenotype.geometry.radius();


		// e.cell_defaults().phenotype.mechanics.attachment_elastic_constant();

		// create a rectangle of cells
		real_t x_offset = -50;
		real_t y_offset = 50;
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
				x_offset = -50;
				y_offset += 2. * element_radius;
			}
		}
	}

	{
		index_t N = e.cell_definitions[1]->phenotype.mechanics.detachment_rate();
		real_t element_radius = e.cell_definitions[1]->phenotype.geometry.radius();


		// e.cell_defaults().phenotype.mechanics.attachment_elastic_constant();

		// create a rectangle of cells
		real_t x_offset = 110;
		real_t y_offset = 0;
		for (index_t i = 0; i < N; i++)
		{
			auto cell = e.get_container().create_cell(*e.cell_definitions[1]);
			cell->residency() = 4;

			point_t<real_t, 3> position = { 0, 0, 0 };

			position[0] = x_offset;
			position[1] = y_offset;
			position[2] = 0;

			cell->assign_position(position);

			std::cout << "placing cell " << i << " at position " << position[0] << " " << position[1] << std::endl;
			x_offset += 2. * element_radius;
			if ((i + 1) % (int)(N / std::sqrt(N)) == 0 && i != 0)
			{
				x_offset = 110;
				y_offset += 2. * element_radius;
			}
		}
	}
}

cell_coloring_funct_t get_my_coloring_function(User_Parameters&) { return paint_by_number_cell_coloring; }
