#include "custom.h"

#include <cmath>

#include "src/original/modules/geometry.h"
#include "src/original/modules/pathology.h"
#include "src/original/modules/settings.h"
#include "src/original/signal_behavior.h"
#include "src/original/standard_models.h"
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
	real_t Xmin = e.m.mesh.bounding_box_mins[0];
	real_t Ymin = e.m.mesh.bounding_box_mins[1];
	real_t Zmin = e.m.mesh.bounding_box_mins[2];

	real_t Xmax = e.m.mesh.bounding_box_maxs[0];
	real_t Ymax = e.m.mesh.bounding_box_maxs[1];
	real_t Zmax = e.m.mesh.bounding_box_maxs[2];

	if (e.m.mesh.dims == 2)
	{
		Zmin = 0.0;
		Zmax = 0.0;
	}

	real_t Xrange = Xmax - Xmin;
	real_t Yrange = Ymax - Ymin;
	real_t Zrange = Zmax - Zmin;

	// create some of each type of cell

	cell* pC;

	for (index_t k = 0; k < e.cell_definitions_count; k++)
	{
		auto& pCD = *e.cell_definitions[k];
		std::cout << "Placing cells of type " << pCD.name << " ... " << std::endl;
		for (index_t n = 0; n < parameters.ints("number_of_cells"); n++)
		{
			point_t<real_t, 3> position = { 0, 0, 0 };
			position[0] = Xmin + random::instance().uniform() * Xrange;
			position[1] = Ymin + random::instance().uniform() * Yrange;
			position[2] = Zmin + random::instance().uniform() * Zrange;

			pC = e.get_container().create_cell(pCD);
			pC->assign_position(position);
		}
	}
	std::cout << std::endl;

	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(config_root, e);
}

cell_coloring_funct_t get_my_coloring_function(User_Parameters&) { return paint_by_number_cell_coloring; }
