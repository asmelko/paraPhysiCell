#include "custom.h"

#include <cmath>

#include "src/original/modules/geometry.h"
#include "src/original/modules/pathology.h"
#include "src/original/modules/settings.h"
#include "src/original/signal_behavior.h"
#include "src/original/standard_models.h"
#include "src/random.h"

std::function<void(cell&)> get_cargo_cell_rule(User_Parameters& parameters)
{
	return [elastic_coefficient = parameters.doubles("elastic_coefficient")](cell& c) {
		set_single_behavior(&c, "cell-cell adhesion elastic constant", elastic_coefficient, c.e());
	};
}

std::function<void(cell&)> get_worker_cell_rule(User_Parameters& parameters)
{
	return [threshold = parameters.doubles("drop_threshold"),
			attached_worker_migration_bias = parameters.doubles("attached_worker_migration_bias"),
			unattached_worker_migration_bias = parameters.doubles("unattached_worker_migration_bias"),
			elastic_coefficient = parameters.doubles("elastic_coefficient")](cell& c) {
		auto* pCell = &c;

		real_t director_signal = get_single_signal(pCell, "director signal", c.e());

		set_single_behavior(pCell, "cell-cell adhesion elastic constant", elastic_coefficient, c.e());

		// have I arrived? If so, release my cargo
		// set chemotaxis weights to seek cargo
		// set migration bias
		if (director_signal > threshold)
		{
			// set receptor = 0 for cells we're detaching from
			// and set their cycle rate to zero
			for (std::size_t k = 0; k < pCell->state.attached_cells().size(); k++)
			{
				cell* pTemp = c.container().get_at(pCell->state.attached_cells()[k]);

				set_single_behavior(pTemp, "custom:receptor", 0.0, c.e());
				set_single_behavior(pTemp, "cycle entry", 0.0, c.e());
			}

			pCell->remove_all_attached_cells();

			set_single_behavior(pCell, "chemotactic response to director signal", 0.0, c.e());
			set_single_behavior(pCell, "chemotactic response to cargo signal", 1.0, c.e());
			set_single_behavior(pCell, "migration bias", unattached_worker_migration_bias, c.e());
		}

		// am I searching for cargo? if so, see if I've found it
		// if( pCell->state.neighbors.size() == 0 ) // pre 1.8.0
		if (pCell->state.attached_cells().size() == 0)
		{
			for (std::size_t i = 0; i < pCell->cells_in_my_mechanics_voxel().size(); i++)
			{
				auto nearby = c.container().get_at(pCell->cells_in_my_mechanics_voxel()[i]);

				// if it is expressing the receptor, dock with it
				// set chemotaxis weights
				// set migration bias

				real_t receptor = get_single_signal(nearby, "custom:receptor", c.e());
				if (receptor > 0.5)
				{
					cell::attach_cells(*pCell, *nearby);
					set_single_behavior(nearby, "custom:receptor", 0.0, c.e());

					set_single_behavior(nearby, "director signal secretion", 0.0, c.e());
					set_single_behavior(nearby, "cargo signal secretion", 0.0, c.e());

					set_single_behavior(pCell, "chemotactic response to director signal", 1.0, c.e());
					set_single_behavior(pCell, "chemotactic response to cargo signal", 0.0, c.e());
					set_single_behavior(pCell, "migration bias", attached_worker_migration_bias, c.e());
				}
			}
		}
	};
}

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

	auto pCD = builder.find_cell_definition("cargo cell");
	pCD->functions.update_phenotype = get_cargo_cell_rule(builder.get_parameters());
	pCD->functions.contact_function = standard_elastic_contact_function;

	pCD = builder.find_cell_definition("worker cell");
	pCD->functions.update_phenotype = get_worker_cell_rule(builder.get_parameters());
	pCD->functions.contact_function = standard_elastic_contact_function;
}

void setup_microenvironment(microenvironment_builder&)
{
	// set domain parameters

	// put any custom code to set non-homogeneous initial conditions or
	// extra Dirichlet nodes here.
}

void create_cargo_cluster_6(const point_t<real_t, 3>& center, environment& e)
{
	// create a hollow cluster at position, with random orientation

	auto pCargoDef = e.find_cell_definition("cargo cell");

	static real_t spacing = 0.95 * pCargoDef->phenotype.geometry.radius() * 2.0;
	static real_t d_Theta = 1.047197551196598; // 2*pi / 6.0

	real_t theta = 6.283185307179586 * random::instance().uniform();

	point_t<real_t, 3> position { 0, 0, 0 };

	cell* pC;
	for (index_t i = 0; i < 6; i++)
	{
		pC = e.cast_container<cell_container>().create_cell(*pCargoDef);

		position[0] = center[0] + spacing * std::cos(theta);
		position[1] = center[1] + spacing * std::sin(theta);

		pC->assign_position(position);

		theta += d_Theta;
	}

	return;
}

void create_cargo_cluster_7(const point_t<real_t, 3>& center, environment& e)
{
	auto pCargoDef = e.find_cell_definition("cargo cell");

	// create a filled cluster at position, with random orientation

	create_cargo_cluster_6(center, e);
	cell* pC = e.cast_container<cell_container>().create_cell(*pCargoDef);
	pC->assign_position(center);

	return;
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

			pC = e.cast_container<cell_container>().create_cell(pCD);
			pC->assign_position(position);
		}
	}
	std::cout << std::endl;

	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(config_root, e);

	/* custom loading here */

	index_t number_of_directors = parameters.ints("number_of_directors");			// 15;
	index_t number_of_cargo_clusters = parameters.ints("number_of_cargo_clusters"); // 100;
	index_t number_of_workers = parameters.ints("number_of_workers");				// 50;

	auto pCargoDef = e.find_cell_definition("cargo cell");
	auto pDirectorDef = e.find_cell_definition("director cell");
	auto pWorkerDef = e.find_cell_definition("worker cell");


	std::cout << "Placing cells ... " << std::endl;

	// randomly place seed cells

	point_t<real_t, 3> position { 0, 0, 0 };

	real_t relative_margin = 0.2;
	real_t relative_outer_margin = 0.02;

	std::cout << "\tPlacing " << number_of_directors << " director cells ... " << std::endl;
	for (index_t i = 0; i < number_of_directors; i++)
	{
		// pick a random location
		position[0] = Xmin + Xrange * (relative_margin + (1.0 - 2 * relative_margin) * random::instance().uniform());

		position[1] =
			Ymin + Yrange * (relative_outer_margin + (1.0 - 2 * relative_outer_margin) * random::instance().uniform());

		// place the cell
		auto pC = e.cast_container<cell_container>().create_cell(*pDirectorDef);
		pC->assign_position(position);
		set_single_behavior(pC, "movable", false, e);
	}

	// place cargo clusters on the fringes

	std::cout << "\tPlacing cargo cells ... " << std::endl;
	for (index_t i = 0; i < number_of_cargo_clusters; i++)
	{
		// pick a random location

		position[0] =
			Xmin + Xrange * (relative_outer_margin + (1 - 2.0 * relative_outer_margin) * random::instance().uniform());

		position[1] =
			Ymin + Yrange * (relative_outer_margin + (1 - 2.0 * relative_outer_margin) * random::instance().uniform());

		if (random::instance().uniform() < 0.5)
		{
			auto pC = e.cast_container<cell_container>().create_cell(*pCargoDef);
			pC->assign_position(position);
		}
		else
		{
			create_cargo_cluster_7(position, e);
		}
	}

	// place "workersworkers"

	std::cout << "\tPlacing worker cells ... " << std::endl;
	for (index_t i = 0; i < number_of_workers; i++)
	{
		// pick a random location

		position[0] = Xmin + Xrange * (relative_margin + (1.0 - 2 * relative_margin) * random::instance().uniform());

		position[1] =
			Ymin + Yrange * (relative_outer_margin + (1.0 - 2 * relative_outer_margin) * random::instance().uniform());

		// place the cell
		auto pC = e.cast_container<cell_container>().create_cell(*pWorkerDef);
		pC->assign_position(position);
	}


	std::cout << "done!" << std::endl;
}

cell_coloring_funct_t get_robot_coloring_function(User_Parameters& parameters)
{
	return [worker_color = parameters.strings("worker_color"), cargo_color = parameters.strings("cargo_color"),
			director_color = parameters.strings("director_color")](cell* pCell) {
		std::string color = "black";
		std::vector<std::string> output(4, color);

		// black cells if necrotic
		if (pCell->phenotype.death.dead() == true)
		{
			return output;
		}

		output[3] = "none"; // no nuclear outline color

		auto worker_ID = pCell->e().find_cell_definition("worker cell")->type;
		auto cargo_ID = pCell->e().find_cell_definition("cargo cell")->type;
		auto director_ID = pCell->e().find_cell_definition("director cell")->type;

		if (pCell->type == worker_ID)
		{
			color = worker_color;
		}
		else if (pCell->type == cargo_ID)
		{
			color = cargo_color;
		}
		else if (pCell->type == director_ID)
		{
			color = director_color;
		}
		else
		{
			color = "white";
		}

		output[0] = color;
		output[2] = color;

		return output;
	};
}
