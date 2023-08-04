#include "custom.h"

#include <cmath>

#include "src/cell_container.h"
#include "src/original/modules/geometry.h"
#include "src/original/modules/settings.h"
#include "src/original/signal_behavior.h"
#include "src/original/standard_models.h"
#include "src/random.h"

void tumor_cell_phenotype_with_oncoprotein(cell& pCell)
{
	update_cell_and_death_parameters_O2_based(pCell);

	// if cell is dead, don't bother with future phenotype changes.
	if (get_single_signal(&pCell, "dead", pCell.e()) > 0.5)
	{
		pCell.functions.update_phenotype = NULL;
		return;
	}

	// multiply proliferation rate by the oncoprotein

	real_t cycle_rate = get_single_behavior(&pCell, "cycle entry", pCell.e());
	cycle_rate *= get_single_signal(&pCell, "custom:oncoprotein", pCell.e());
	set_single_behavior(&pCell, "cycle entry", cycle_rate, pCell.e());
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

	auto pCD = builder.find_cell_definition("cancer cell");
	pCD->functions.update_phenotype = tumor_cell_phenotype_with_oncoprotein;

	pCD->parameters.o2_proliferation_saturation = 38;
	pCD->parameters.o2_reference = 38;
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
		auto pCD = e.cell_definitions[k];
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

	// custom placement

	auto pCD = e.find_cell_definition("cancer cell");
	real_t cell_radius = pCD->phenotype.geometry.radius();
	real_t cell_spacing = 0.95 * 2.0 * cell_radius;

	real_t tumor_radius = parameters.doubles("tumor_radius"); // 250.0;

	// Parameter<real_t> temp;

	// index_t i = parameters.doubles.find_index("tumor_radius");

	cell* pCell = NULL;

	real_t x = 0.0;
	real_t x_outer = tumor_radius;
	real_t y = 0.0;

	real_t p_mean = parameters.doubles("oncoprotein_mean");
	real_t p_sd = parameters.doubles("oncoprotein_sd");
	real_t p_min = parameters.doubles("oncoprotein_min");
	real_t p_max = parameters.doubles("oncoprotein_max");

	auto& container = e.cast_container<cell_container>();

	index_t n = 0;
	while (y < tumor_radius)
	{
		x = 0.0;
		if (n % 2 == 1)
		{
			x = 0.5 * cell_spacing;
		}
		x_outer = sqrt(tumor_radius * tumor_radius - y * y);

		while (x < x_outer)
		{
			pCell = container.create_cell(*pCD); // tumor cell
			pCell->assign_position({ x, y, 0.0 });
			real_t p = random::instance().normal(p_mean, p_sd);
			if (p < p_min)
			{
				p = p_min;
			}
			if (p > p_max)
			{
				p = p_max;
			}
			set_single_behavior(pCell, "custom:oncoprotein", p, e);

			if (fabs(y) > 0.01)
			{
				pCell = container.create_cell(*pCD); // tumor cell
				pCell->assign_position({ x, -y, 0.0 });
				real_t p = random::instance().normal(p_mean, p_sd);
				if (p < p_min)
				{
					p = p_min;
				}
				if (p > p_max)
				{
					p = p_max;
				}
				set_single_behavior(pCell, "custom:oncoprotein", p, e);
			}

			if (fabs(x) > 0.01)
			{
				pCell = container.create_cell(*pCD); // tumor cell
				pCell->assign_position({ -x, y, 0.0 });
				real_t p = random::instance().normal(p_mean, p_sd);
				if (p < p_min)
				{
					p = p_min;
				}
				if (p > p_max)
				{
					p = p_max;
				}
				set_single_behavior(pCell, "custom:oncoprotein", p, e);

				if (fabs(y) > 0.01)
				{
					pCell = container.create_cell(*pCD); // tumor cell
					pCell->assign_position({ -x, -y, 0.0 });
					real_t p = random::instance().normal(p_mean, p_sd);
					if (p < p_min)
					{
						p = p_min;
					}
					if (p > p_max)
					{
						p = p_max;
					}
					set_single_behavior(pCell, "custom:oncoprotein", p, e);
				}
			}
			x += cell_spacing;
		}

		y += cell_spacing * sqrt(3.0) / 2.0;
		n++;
	}

	real_t sum = 0.0;
	real_t min = 9e9;
	real_t max = -9e9;
	for (index_t i = 0; i < container.agents_count(); i++)
	{
		real_t r = get_single_signal(container.get_at(i), "custom:oncoprotein", e);
		sum += r;
		if (r < min)
		{
			min = r;
		}
		if (r > max)
		{
			max = r;
		}
	}
	real_t mean = sum / (container.agents_count() + 1e-15);
	// compute standard deviation
	sum = 0.0;
	for (index_t i = 0; i < container.agents_count(); i++)
	{
		real_t r = get_single_signal(container.get_at(i), "custom:oncoprotein", e);
		sum += (r - mean) * (r - mean);
	}
	real_t standard_deviation = sqrt(sum / (container.agents_count() - 1.0 + 1e-15));

	std::cout << std::endl << "Oncoprotein summary: " << std::endl << "===================" << std::endl;
	std::cout << "mean: " << mean << std::endl;
	std::cout << "standard deviation: " << standard_deviation << std::endl;
	std::cout << "[min max]: [" << min << " " << max << "]" << std::endl << std::endl;

	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(config_root, e);

	return;
}

cell_coloring_funct_t get_heterogeneity_coloring_function(User_Parameters& parameters)
{
	return [p_min = parameters.doubles("oncoprotein_min"), p_max = parameters.doubles("oncoprotein_max")](cell* pCell) {
		real_t p = get_single_signal(pCell, "custom:oncoprotein", pCell->e());

		// immune are black
		std::vector<std::string> output(4, "black");

		if (pCell->type == 1)
		{
			return output;
		}

		// live cells are green, but shaded by oncoprotein value
		if (pCell->phenotype.death.dead() == false)
		{
			index_t oncoprotein = (index_t)round((1.0 / (p_max - p_min)) * (p - p_min) * 255.0);
			char szTempString[128];
			sprintf(szTempString, "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255 - oncoprotein);
			output[0].assign(szTempString);
			output[1].assign(szTempString);

			sprintf(szTempString, "rgb(%u,%u,%u)", (index_t)round(output[0][0] / p_max),
					(index_t)round(output[0][1] / p_max), (index_t)round(output[0][2] / p_max));
			output[2].assign(szTempString);

			return output;
		}

		// if not, dead colors

		if (get_single_signal(pCell, "apoptotic", pCell->e()) > 0.5)
		{
			output[0] = "rgb(255,0,0)";
			output[2] = "rgb(125,0,0)";
		}

		// Necrotic - Brown
		if (get_single_signal(pCell, "necrotic", pCell->e()) > 0.5)
		{
			output[0] = "rgb(250,138,38)";
			output[2] = "rgb(139,69,19)";
		}

		return output;
	};
}
