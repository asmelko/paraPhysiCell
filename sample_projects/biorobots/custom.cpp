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

		double director_signal = get_single_signal(pCell, "director signal", c.e());
		double cargo_signal = get_single_signal(pCell, "cargo signal", c.e());

		set_single_behavior(pCell, "cell-cell adhesion elastic constant", elastic_coefficient, c.e());

		// have I arrived? If so, release my cargo
		// set chemotaxis weights to seek cargo
		// set migration bias
		if (director_signal > threshold)
		{
			// set receptor = 0 for cells we're detaching from
			// and set their cycle rate to zero
			for (int k = 0; k < pCell->state.attached_cells().size(); k++)
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
			for (int i = 0; i < pCell->state.neighbors().size(); i++)
			{
				auto nearby = c.container().get_at(pCell->state.neighbors()[i]);

				// if it is expressing the receptor, dock with it
				// set chemotaxis weights
				// set migration bias

				double receptor = get_single_signal(nearby, "custom:receptor", c.e());
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
	auto& cell_definitions = builder.get_cell_definitions();

	/*
	   Put any modifications to individual cell definitions here.

	   This is a good place to set custom functions.
	*/

	auto pCD = builder.find_cell_definition("cargo cell");
	pCD->functions.update_phenotype = get_cargo_cell_rule(builder.get_parameters());
	pCD->functions.contact_function = standard_elastic_contract_function;

	pCD = builder.find_cell_definition("worker cell");
	pCD->functions.update_phenotype = get_worker_cell_rule(builder.get_parameters());
	pCD->functions.contact_function = standard_elastic_contract_function;
}

void setup_microenvironment(microenvironment_builder& m_builder)
{
	// set domain parameters

	// put any custom code to set non-homogeneous initial conditions or
	// extra Dirichlet nodes here.
}

void create_cargo_cluster_6(const point_t<real_t, 3>& center, environment& e)
{
	// create a hollow cluster at position, with random orientation

	auto pCargoDef = e.find_cell_definition("cargo cell");

	static double spacing = 0.95 * pCargoDef->phenotype.geometry.radius() * 2.0;
	static double d_Theta = 1.047197551196598; // 2*pi / 6.0

	double theta = 6.283185307179586 * random::instance().uniform();

	point_t<real_t, 3> position { 0, 0, 0 };

	cell* pC;
	for (int i = 0; i < 6; i++)
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
	double Xmin = e.m.mesh.bounding_box_mins[0];
	double Ymin = e.m.mesh.bounding_box_mins[1];
	double Zmin = e.m.mesh.bounding_box_mins[2];

	double Xmax = e.m.mesh.bounding_box_maxs[0];
	double Ymax = e.m.mesh.bounding_box_maxs[1];
	double Zmax = e.m.mesh.bounding_box_maxs[2];

	if (e.m.mesh.dims == 2)
	{
		Zmin = 0.0;
		Zmax = 0.0;
	}

	double Xrange = Xmax - Xmin;
	double Yrange = Ymax - Ymin;
	double Zrange = Zmax - Zmin;

	// create some of each type of cell

	cell* pC;

	for (int k = 0; k < e.cell_definitions_count; k++)
	{
		auto pCD = e.cell_definitions[k];
		std::cout << "Placing cells of type " << pCD.name << " ... " << std::endl;
		for (int n = 0; n < parameters.ints("number_of_cells"); n++)
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

	int number_of_directors = parameters.ints("number_of_directors");			// 15;
	int number_of_cargo_clusters = parameters.ints("number_of_cargo_clusters"); // 100;
	int number_of_workers = parameters.ints("number_of_workers");				// 50;

	auto pCargoDef = e.find_cell_definition("cargo cell");
	auto pDirectorDef = e.find_cell_definition("director cell");
	auto pWorkerDef = e.find_cell_definition("worker cell");


	std::cout << "Placing cells ... " << std::endl;

	// randomly place seed cells

	point_t<real_t, 3> position { 0, 0, 0 };

	double relative_margin = 0.2;
	double relative_outer_margin = 0.02;

	std::cout << "\tPlacing " << number_of_directors << " director cells ... " << std::endl;
	for (int i = 0; i < number_of_directors; i++)
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
	for (int i = 0; i < number_of_cargo_clusters; i++)
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
	for (int i = 0; i < number_of_workers; i++)
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

std::vector<std::string> my_coloring_function(cell* pCell) { return paint_by_number_cell_coloring(pCell); }
