#include "custom.h"

#include "src/original/core/basic_signaling.h"
#include "src/original/core/constants.h"
#include "src/original/modules/geometry.h"
#include "src/original/modules/pathology.h"
#include "src/random.h"

void bacteria_phenotype(cell& pCell)
{
	// find my cell definition
	auto pCD = pCell.e().find_cell_definition(pCell.type_name);

	// sample resource, quorum, and toxin

	auto nR = pCell.e().m.find_substrate_index("resource");
	auto nDebris = pCell.e().m.find_substrate_index("debris");
	auto nQuorum = pCell.e().m.find_substrate_index("quorum");
	auto nToxin = pCell.e().m.find_substrate_index("toxin");

	// if dead: stop exporting quorum factor.
	// also, replace phenotype function
	if (pCell.phenotype.death.dead() == true)
	{
		pCell.phenotype.secretion.net_export_rates()[nQuorum] = 0;
		pCell.phenotype.secretion.net_export_rates()[nToxin] = 0;

		pCell.phenotype.secretion.net_export_rates()[nDebris] = pCell.phenotype.volume.total();

		pCell.functions.update_phenotype = NULL;

		return;
	}

	auto samples = pCell.nearest_density_vector();
	real_t R = samples[nR];
	real_t Q = samples[nQuorum];

	// resource increases cycle entry
	real_t base_val = pCD->phenotype.cycle.data.exit_rate(0);
	real_t max_val = base_val * 10.0;
	real_t min_cycle_resource = pCD->custom_data["cycling_entry_threshold_resource"]; // 0.15
	pCell.phenotype.cycle.data.exit_rate(0) = max_val * linear_response_function(R, min_cycle_resource, 1);

	// resource decreses necrosis

	max_val = 0.0028;
	index_t nNecrosis = pCell.phenotype.death.find_death_model_index(constants::necrosis_death_model);
	real_t saturation_necrosis_resource = pCD->custom_data["necrosis_saturation_resource"]; // 0.075
	real_t threshold_necrosis_resource = pCD->custom_data["necrosis_threshold_resource"];	// 0.15
	pCell.phenotype.death.rates[nNecrosis] =
		max_val * decreasing_linear_response_function(R, saturation_necrosis_resource, threshold_necrosis_resource);

	// resource decreases motile speed

	real_t signal = R;
	base_val = pCD->phenotype.motility.migration_speed();
	real_t max_response = 0.0;
	real_t motility_resource_halfmax =
		pCD->custom_data["migration_speed_halfmax"]; // 0.25 //
													 // parameters.doubles("bacteria_motility_resource_halfmax");
	real_t hill = Hill_response_function(signal, motility_resource_halfmax, 1.5);
	pCell.phenotype.motility.migration_speed() = base_val + (max_response - base_val) * hill;

	// quorum and resource increases motility bias
	signal = Q + R;
	base_val = pCD->phenotype.motility.migration_speed();
	max_response = 1.0;
	real_t bias_halfmax = pCD->custom_data["migration_bias_halfmax"];
	// 0.5 //  parameters.doubles("bacteria_migration_bias_halfmax");
	hill = Hill_response_function(signal, bias_halfmax, 1.5);
	pCell.phenotype.motility.migration_bias() = base_val + (max_response - base_val) * hill;

	// damage increases death
	index_t nApoptosis = pCell.phenotype.death.find_death_model_index(constants::apoptosis_death_model);

	signal = pCell.phenotype.cell_integrity.damage();
	base_val = pCD->phenotype.death.rates[nApoptosis];

	real_t damage_halfmax = pCD->custom_data["damage_halfmax"];
	real_t relative_max_damage_death = pCD->custom_data["relative_max_damage_death"];
	max_response = base_val * relative_max_damage_death;

	// 36 // parameters.doubles("bacteria_damage_halfmax");
	hill = Hill_response_function(signal, damage_halfmax, 1.5);
	pCell.phenotype.death.rates[nApoptosis] = base_val + (max_response - base_val) * hill;
}

/* https://www.karger.com/Article/Fulltext/494069 */

void macrophage_phenotype(cell& pCell)
{
	// find my cell definition
	auto pCD = pCell.e().find_cell_definition(pCell.type_name);

	// sample environment

	index_t nPIF = pCell.e().m.find_substrate_index("pro-inflammatory");
	index_t nDebris = pCell.e().m.find_substrate_index("debris");
	index_t nQ = pCell.e().m.find_substrate_index("quorum");

	// if dead, release debris
	if (pCell.phenotype.death.dead() == true)
	{
		pCell.phenotype.secretion.net_export_rates()[nDebris] = pCell.phenotype.volume.total();
		pCell.functions.update_phenotype = NULL;
		return;
	}

	auto samples = pCell.nearest_density_vector();
	real_t debris = samples[nDebris];
	real_t Q = samples[nQ];

	// sample contacts

	auto bacteria_type = pCell.e().find_cell_definition("bacteria")->type;

	index_t num_bacteria = 0;
	index_t num_dead = 0;
	for (std::size_t n = 0; n < pCell.state.neighbors().size(); n++)
	{
		cell* pC = pCell.container().get_at(pCell.state.neighbors()[n]);
		if (pC->phenotype.death.dead() == true)
		{
			num_dead++;
		}
		else
		{
			if (pC->type == bacteria_type)
			{
				num_bacteria++;
			}
		}
	}

	// contact with dead cells or bacteria, or debris
	// increases secretion of pro-inflammatory

	real_t secretion_dead_sensitivity = 1;
	real_t secretion_bacteria_sensitivity = 1;
	real_t secretion_debris_sensitivity = 2;
	real_t secretion_quorum_sensitivity = 5;

	real_t base_val = pCD->phenotype.secretion.secretion_rates()[nPIF];
	real_t max_response = 10; // pCell.phenotype.volume.total();
	real_t signal = secretion_dead_sensitivity * num_dead + secretion_bacteria_sensitivity * num_bacteria
					+ secretion_debris_sensitivity * debris + secretion_quorum_sensitivity * Q;
	real_t half_max = pCD->custom_data["secretion_halfmax"]; // 0.5; // 0.5;
	real_t hill = Hill_response_function(signal, half_max, 1.5);


	pCell.phenotype.secretion.secretion_rates()[nPIF] = base_val + (max_response - base_val) * hill;

	/*
		#pragma omp critical
		{
		std::cout << "secretion index: " << nPIF << " base: " << base_val << " max: " << max_response << " actual: " <<
	   phenotype.secretion.secretion_rates[nPIF] << std::endl; std::cout << "\tsignal: " << signal << " vs halfmax: " <<
	   half_max << std::endl; std::cout << "\t\tdead: " << num_dead << " bac: " << num_bacteria << " debris: " << debris
	   << " Q: " << Q << std::endl; std::cout << "\t\t\tsaturation: " <<
	   phenotype.secretion.saturation_densities[nPIF]<< std::endl;
		}
	*/

	// chemotaxis bias increases with debris or quorum factor

	real_t bias_debris_sensitivity = 0.1;
	real_t bias_quorum_sensitivity = 1;

	base_val = pCD->phenotype.motility.migration_bias();
	max_response = 0.75;
	signal = bias_debris_sensitivity * debris + bias_quorum_sensitivity * Q; // + 10 * PIF;
	half_max = pCD->custom_data["migration_bias_halfmax"];					 // 0.01 // 0.005 //0.1 // 0.05
	hill = Hill_response_function(signal, half_max, 1.5);
	pCell.phenotype.motility.migration_bias() = base_val + (max_response - base_val) * hill;

	/*
		#pragma omp critical
		{
		std::cout << "signal: " << signal << " halfmax: " << half_max
		<< " hill: " << hill << std::endl;

		std::cout << "\tbase: " << base_val
		<< " max: " << max_response
		<< " actual: " << phenotype.motility.migration_bias << std::endl;
		}
	*/

	// migration speed slows down in the presence of debris or quorum factor

	base_val = pCD->phenotype.motility.migration_speed();
	max_response = 0.1 * base_val;
	signal = bias_debris_sensitivity * debris + bias_quorum_sensitivity * Q; // + 10 * PIF;
	half_max = pCD->custom_data["migration_speed_halfmax"];					 // 0.1 // 0.05
	hill = Hill_response_function(signal, half_max, 1.5);
	pCell.phenotype.motility.migration_speed() = base_val + (max_response - base_val) * hill;

	return;
}

void CD8Tcell_phenotype(cell& pCell)
{
	// find my cell definition
	auto pCD = pCell.e().find_cell_definition(pCell.type_name);

	// sample environment

	index_t nDebris = pCell.e().m.find_substrate_index("debris");
	index_t nPIF = pCell.e().m.find_substrate_index("pro-inflammatory");

	auto samples = pCell.nearest_density_vector();
	real_t PIF = samples[nPIF];

	// if dead, release debris
	if (pCell.phenotype.death.dead() == true)
	{
		pCell.phenotype.secretion.net_export_rates()[nDebris] = pCell.phenotype.volume.total();
		pCell.functions.update_phenotype = NULL;
		return;
	}

	// migration bias increases with pro-inflammatory

	real_t base_val = pCD->phenotype.motility.migration_bias();
	real_t max_val = 0.75;
	real_t half_max = pCD->custom_data["migration_bias_halfmax"]; // 0.05 // 0.25
	real_t hill = Hill_response_function(PIF, half_max, 1.5);

	pCell.phenotype.motility.migration_bias() = base_val + (max_val - base_val) * hill;

	/*
		#pragma omp critical
		{
			std::cout << "signal: " << signal << " halfmax: " << half_max
			<< " hill: " << hill << std::endl;

			std::cout << "\tbase: " << base_val
			<< " max: " << max_val
			<< " actual: " << phenotype.motility.migration_bias << std::endl;
		}
	*/

	return;
}

void neutrophil_phenotype(cell& pCell)
{
	// find my cell definition
	auto pCD = pCell.e().find_cell_definition(pCell.type_name);

	// sample environment

	index_t nDebris = pCell.e().m.find_substrate_index("debris");
	index_t nPIF = pCell.e().m.find_substrate_index("pro-inflammatory");

	auto samples = pCell.nearest_density_vector();
	real_t PIF = samples[nPIF];

	// if dead, release debris
	if (pCell.phenotype.death.dead() == true)
	{
		pCell.phenotype.secretion.net_export_rates()[nDebris] = pCell.phenotype.volume.total();
		pCell.functions.update_phenotype = NULL;
		return;
	}

	// migration bias increases with pro-inflammatory

	real_t base_val = pCD->phenotype.motility.migration_bias();
	real_t max_val = 0.75;
	real_t half_max = pCD->custom_data["migration_bias_halfmax"]; // 0.25
	real_t hill = Hill_response_function(PIF, half_max, 1.5);

	pCell.phenotype.motility.migration_bias() = base_val + (max_val - base_val) * hill;

	return;
}

std::function<void(cell&)> get_stem_cell_phenotype(User_Parameters& parameters)
{
	return [max_stem_diff = parameters.doubles("max_stem_differentiation")](cell& pCell) {
		// find my cell definition
		auto pCD = pCell.e().find_cell_definition(pCell.type_name);

		// sample environment

		index_t nR = pCell.e().m.find_substrate_index("resource");
		index_t nTox = pCell.e().m.find_substrate_index("toxin");
		index_t nDebris = pCell.e().m.find_substrate_index("debris");

		// if dead, release debris
		if (pCell.phenotype.death.dead() == true)
		{
			pCell.phenotype.secretion.net_export_rates()[nDebris] = pCell.phenotype.volume.total();
			pCell.functions.update_phenotype = NULL;
			return;
		}

		auto samples = pCell.nearest_density_vector();
		real_t R = samples[nR];
		real_t toxin = samples[nTox];

		// sample contacts

		index_t stem_type = pCell.e().find_cell_definition("stem")->type;
		index_t diff_type = pCell.e().find_cell_definition("differentiated")->type;

		index_t num_stem = 0;
		index_t num_differentiated = 0;
		for (std::size_t n = 0; n < pCell.state.neighbors().size(); n++)
		{
			cell* pC = pCell.container().get_at(pCell.state.neighbors()[n]);
			if (pC->type == stem_type)
			{
				num_stem++;
			}
			if (pC->type == num_differentiated)
			{
				num_differentiated++;
			}
		}

		// contact with a stem cell increases differentiation
		real_t stem_diff_halfmax = pCD->custom_data["differentiation_contact_halfmax"]; // 0.1

		real_t base_val = 0;			// phenotype.cell_transformations.transformation_rates[diff_type];
		real_t max_val = max_stem_diff; // 0.0075;
		real_t signal = num_stem;
		real_t half_max = stem_diff_halfmax; // 0.1;
		real_t hill = Hill_response_function(signal, half_max, 1.5);
		pCell.phenotype.cell_transformations.transformation_rates()[diff_type] = base_val + (max_val - base_val) * hill;

		// contact with a differentiated cell reduces proliferation
		// high rate of proliferation unless in contact with a differentiated cell

		real_t stem_cycling_halfmax = pCD->custom_data["cycling_contact_halfmax"]; // 0.1;

		base_val = pCD->phenotype.cycle.data.exit_rate(0); // 0.002;
		max_val = 0.0;
		signal = num_differentiated;
		half_max = stem_cycling_halfmax; //  0.1;
		hill = Hill_response_function(signal, half_max, 1.5);
		pCell.phenotype.cycle.data.exit_rate(0) = base_val + (max_val - base_val) * hill;

		// resource reduces necrotic death

		max_val = 0.0028;
		index_t nNecrosis = pCell.phenotype.death.find_death_model_index(constants::necrosis_death_model);
		real_t stem_saturation_necrosis = pCD->custom_data["necrosis_saturation_resource"];
		real_t stem_threshold_necrosis = pCD->custom_data["necrosis_threshold_resource"];

		pCell.phenotype.death.rates[nNecrosis] =
			max_val * decreasing_linear_response_function(R, stem_saturation_necrosis, stem_threshold_necrosis);

		// toxin increases apoptotic death

		index_t nApoptosis = pCell.phenotype.death.find_death_model_index(constants::apoptosis_death_model);

		real_t toxicity_halfmax = pCD->custom_data["toxicity_halfmax"]; // 0.4
		real_t relative_max_toxicity = pCD->custom_data["relative_max_toxicity"];

		signal = toxin;
		base_val = pCD->phenotype.death.rates[nApoptosis];
		max_val = base_val * relative_max_toxicity; // 100*base_val;

		hill = Hill_response_function(signal, toxicity_halfmax, 1.5);
		pCell.phenotype.death.rates[nApoptosis] = base_val + (max_val - base_val) * hill;
	};
}

void differentiated_cell_phenotype(cell& pCell)
{
	// find my cell definition
	auto pCD = pCell.e().find_cell_definition(pCell.type_name);

	// sample environment

	index_t nR = pCell.e().m.find_substrate_index("resource");
	index_t nTox = pCell.e().m.find_substrate_index("toxin");
	index_t nDebris = pCell.e().m.find_substrate_index("debris");

	// if dead, release debris
	if (pCell.phenotype.death.dead() == true)
	{
		pCell.phenotype.secretion.net_export_rates()[nDebris] = pCell.phenotype.volume.total();
		pCell.functions.update_phenotype = NULL;
		return;
	}

	auto samples = pCell.nearest_density_vector();
	real_t R = samples[nR];
	real_t toxin = samples[nTox];


	real_t signal = 0.0;
	real_t hill = 0.0;

	// pressure reduces proliferation
	signal = pCell.state.simple_pressure();
	real_t pressure_halfmax = pCD->custom_data["cycling_pressure_halfmax"]; // 0.5
	hill = Hill_response_function(signal, pressure_halfmax, 1.5);
	real_t base_val = pCD->phenotype.cycle.data.exit_rate(0);

	pCell.phenotype.cycle.data.exit_rate(0) = (1 - hill) * base_val;

	// resource reduces necrotic death

	real_t max_val = 0.0028;
	index_t nNecrosis = pCell.phenotype.death.find_death_model_index(constants::necrosis_death_model);

	// get same from bacteria
	real_t necrosis_saturation = pCD->custom_data["necrosis_saturation_resource"]; // 0.075
	real_t necrosis_threshold = pCD->custom_data["necrosis_threshold_resource"];   // 0.15

	pCell.phenotype.death.rates[nNecrosis] =
		max_val * decreasing_linear_response_function(R, necrosis_saturation, necrosis_threshold);

	// toxin increases apoptotic death

	index_t nApoptosis = pCell.phenotype.death.find_death_model_index(constants::apoptosis_death_model);

	real_t toxicity_halfmax = pCD->custom_data["toxicity_halfmax"];			   // 0.2
	real_t relative_max_tox_death = pCD->custom_data["relative_max_toxicity"]; // 100

	signal = toxin;
	base_val = pCD->phenotype.death.rates[nApoptosis];
	real_t max_response = base_val * relative_max_tox_death;
	hill = Hill_response_function(signal, toxicity_halfmax, 1.5);
	// std::cout << "tox: " << signal << " " << hill << std::endl;
	pCell.phenotype.death.rates[nApoptosis] = base_val + (max_response - base_val) * hill;

	return;
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

	// set up bacteria

	auto pCD = builder.find_cell_definition("bacteria");
	pCD->functions.update_phenotype = bacteria_phenotype;
	// pCD->functions.update_migration_bias = advanced_chemotaxis_function;
	// pCD->phenotype.motility.chemotactic_sensitivity( "resource" ) = 1;
	// pCD->phenotype.motility.chemotactic_sensitivity( "quorum" ) = 0.1;

	// set up blood vessels

	pCD = builder.find_cell_definition("blood vessel");
	pCD->is_movable = false;

	// set up stem cells

	pCD = builder.find_cell_definition("stem");
	pCD->functions.update_phenotype = get_stem_cell_phenotype(builder.get_parameters());
	// pCD->phenotype.cell_transformations.transformation_rate("differentiated") = 0.0001;

	// set up differentiated cells

	pCD = builder.find_cell_definition("differentiated");
	pCD->functions.update_phenotype = differentiated_cell_phenotype;

	// set up macrophages

	pCD = builder.find_cell_definition("macrophage");
	// pCD->phenotype.cell_index_teractions.dead_phagocytosis_rate = 0.05;
	pCD->functions.update_phenotype = macrophage_phenotype;
	// pCD->functions.update_migration_bias = advanced_chemotaxis_function;
	// pCD->phenotype.motility.chemotactic_sensitivity( "debris" ) = 0.1;
	// pCD->phenotype.motility.chemotactic_sensitivity( "quorum" ) = 1;


	// set up CD8+ T cells
	pCD = builder.find_cell_definition("CD8+ T cell");
	pCD->functions.update_phenotype = CD8Tcell_phenotype;
	// pCD->phenotype.cell_index_teractions.attack_rate("bacteria") = 0.05;

	// set up neutrophil
	pCD = builder.find_cell_definition("neutrophil");
	pCD->functions.update_phenotype = neutrophil_phenotype;
	// pCD->phenotype.cell_index_teractions.live_phagocytosis_rate("bacteria") = 0.05;
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
		index_t num_cells = parameters.ints("number_of_cells");
		if (num_cells > 0)
		{
			std::cout << "Placing cells of type " << pCD.name << " ... " << std::endl;
		}
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

	// parameter-based placement
	// bacteria
	auto pCD = e.find_cell_definition("bacteria");
	std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
	for (index_t n = 0; n < parameters.ints("number_of_bacteria"); n++)
	{
		point_t<real_t, 3> position = { 0, 0, 0 };
		position[0] = Xmin + random::instance().uniform() * Xrange;
		position[1] = Ymin + random::instance().uniform() * Yrange;
		position[2] = Zmin + random::instance().uniform() * Zrange;

		pC = e.get_container().create_cell(*pCD);
		pC->assign_position(position);
	}

	// blood vessels
	pCD = e.find_cell_definition("blood vessel");
	std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
	for (index_t n = 0; n < parameters.ints("number_of_blood_vessels"); n++)
	{
		point_t<real_t, 3> position = { 0, 0, 0 };
		position[0] = Xmin + random::instance().uniform() * Xrange;
		position[1] = Ymin + random::instance().uniform() * Yrange;
		position[2] = Zmin + random::instance().uniform() * Zrange;

		pC = e.get_container().create_cell(*pCD);
		pC->assign_position(position);
	}

	// stem cells
	pCD = e.find_cell_definition("stem");
	std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
	for (index_t n = 0; n < parameters.ints("number_of_stem_cells"); n++)
	{
		point_t<real_t, 3> position = { 0, 0, 0 };
		position[0] = Xmin + random::instance().uniform() * Xrange;
		position[1] = Ymin + random::instance().uniform() * Yrange;
		position[2] = Zmin + random::instance().uniform() * Zrange;

		pC = e.get_container().create_cell(*pCD);
		pC->assign_position(position);
	}

	// differentiated cells
	pCD = e.find_cell_definition("differentiated");
	std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
	for (index_t n = 0; n < parameters.ints("number_of_differentiated_cells"); n++)
	{
		point_t<real_t, 3> position = { 0, 0, 0 };
		position[0] = Xmin + random::instance().uniform() * Xrange;
		position[1] = Ymin + random::instance().uniform() * Yrange;
		position[2] = Zmin + random::instance().uniform() * Zrange;

		pC = e.get_container().create_cell(*pCD);
		pC->assign_position(position);
	}

	// macrophages
	pCD = e.find_cell_definition("macrophage");
	std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
	for (index_t n = 0; n < parameters.ints("number_of_macrophages"); n++)
	{
		point_t<real_t, 3> position = { 0, 0, 0 };
		position[0] = Xmin + random::instance().uniform() * Xrange;
		position[1] = Ymin + random::instance().uniform() * Yrange;
		position[2] = Zmin + random::instance().uniform() * Zrange;

		pC = e.get_container().create_cell(*pCD);
		pC->assign_position(position);
	}

	// neutrophils
	pCD = e.find_cell_definition("neutrophil");
	std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
	for (index_t n = 0; n < parameters.ints("number_of_neutrophils"); n++)
	{
		point_t<real_t, 3> position = { 0, 0, 0 };
		position[0] = Xmin + random::instance().uniform() * Xrange;
		position[1] = Ymin + random::instance().uniform() * Yrange;
		position[2] = Zmin + random::instance().uniform() * Zrange;

		pC = e.get_container().create_cell(*pCD);
		pC->assign_position(position);
	}

	// CD8+ T cells
	pCD = e.find_cell_definition("CD8+ T cell");
	std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
	for (index_t n = 0; n < parameters.ints("number_of_CD8T_cells"); n++)
	{
		point_t<real_t, 3> position = { 0, 0, 0 };
		position[0] = Xmin + random::instance().uniform() * Xrange;
		position[1] = Ymin + random::instance().uniform() * Yrange;
		position[2] = Zmin + random::instance().uniform() * Zrange;

		pC = e.get_container().create_cell(*pCD);
		pC->assign_position(position);
	}

	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(config_root, e);

	return;
}

std::vector<std::string> my_coloring_function(cell* pCell) { return paint_by_number_cell_coloring(pCell); }

std::vector<std::string> my_coloring_function_for_substrate(real_t concentration, real_t max_conc, real_t min_conc)
{
	return paint_by_density_percentage(concentration, max_conc, min_conc);
}
