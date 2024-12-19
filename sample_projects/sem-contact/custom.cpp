#include "custom.h"

#include <cmath>

#include "../helpers.h"
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

	setup_potential_parameters(e, parameters);

	e.cell_definitions[1]->functions.custom_cell_rule = [](cell& c) {
		// forced to move in a direction
		c.velocity()[1] -= 1;
	};

	make_circle(e, e.cell_definitions[2].get(), 0, { -50, 100 });
	make_circle(e, e.cell_definitions[2].get(), 1, { 100, 100 });
	make_packed_square(e, e.cell_definitions[0].get(), 2, { 300, 100 });
	make_packed_square(e, e.cell_definitions[0].get(), 3, { 450, 100 });

	make_circle(e, e.cell_definitions[1].get(), 4, { -50, 200 });
	make_circle(e, e.cell_definitions[1].get(), 5, { 120, 200 });
	make_packed_square(e, e.cell_definitions[1].get(), 6, { 300, 200 });
	make_packed_square(e, e.cell_definitions[1].get(), 7, { 470, 200 });
}

cell_coloring_funct_t get_my_coloring_function(User_Parameters&) { return paint_by_number_cell_coloring; }
