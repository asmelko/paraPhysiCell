#include "builder.h"

#include <iostream>
#include <stdexcept>

#include "../original/constants.h"
#include "../original/modules/pugixml_helper.h"
#include "../original/modules/vector_utils.h"
#include "../original/rules.h"
#include "../original/signal_behavior.h"
#include "../original/standard_models.h"
#include "../random.h"
#include "microenvironment.h"
#include "microenvironment_builder.h"
#include "types.h"

using namespace biofvm;
using namespace physicell;

builder::builder(const std::string& config_path) : config_path_(config_path) { get_config_root(); }

const pugi::xml_node& builder::get_config_root()
{
	std::cout << "Using config file " << config_path_ << " ... " << std::endl;

	pugi::xml_document config_doc;
	pugi::xml_parse_result result = config_doc.load_file(config_path_.c_str());

	if (result.status != pugi::xml_parse_status::status_ok)
	{
		throw std::runtime_error("Error loading " + config_path_);
	}

	config_root_ = config_doc.child("PhysiCell_settings");

	return config_root_;
}

PhysiCell_Settings& builder::get_settings()
{
	if (!settings_)
		settings_->read_from_pugixml(config_root_);
	return *settings_;
}

User_Parameters& builder::get_parameters()
{
	if (!parameters_)
		parameters_->read_from_pugixml(config_root_);
	return *parameters_;
}

microenvironment_builder& builder::get_microenvironment_builder()
{
	auto node = xml_find_node(config_root_, "domain");

	index_t xmin = xml_get_int_value(node, "x_min");
	index_t xmax = xml_get_int_value(node, "x_max");
	index_t ymin = xml_get_int_value(node, "y_min");
	index_t ymax = xml_get_int_value(node, "y_max");
	index_t zmin = xml_get_int_value(node, "z_min");
	index_t zmax = xml_get_int_value(node, "z_max");
	index_t dx = xml_get_int_value(node, "dx");
	index_t dy = xml_get_int_value(node, "dy");
	index_t dz = xml_get_int_value(node, "dz");

	bool simulate_2D = xml_get_bool_value(node, "use_2D");

	index_t dims = 3;

	if (simulate_2D == true)
	{
		zmin = -dz / 2;
		zmax = dz / 2;

		dims = 2;
	}

	m_builder_.resize(dims, { xmin, ymin, zmin }, { xmax, ymax, zmax }, { dx, dy, dz });

	{
		auto& s = get_settings();

		m_builder_.set_space_units(s.space_units);
		m_builder_.set_time_units(s.time_units);

		auto node = xml_find_node(config_root_, "overall");

		// check to see if dt is specified in overall options
		// if so, set from XML

		pugi::xml_node search_result;
		search_result = xml_find_node(node, "dt_diffusion");
		if (search_result)
		{
			m_builder_.set_time_step(xml_get_my_double_value(search_result));
		}
	}

	// First, look for the correct XML node.
	// If it isn't there, return false.

	node = xml_find_node(config_root_, "microenvironment_setup");
	if (!node)
	{
		std::cout << "Warning: microenvironment_setup not found in config file." << std::endl;
		return m_builder_;
	}

	// now that we're using the XML to specify the microenvironment, don't
	// use old defaults

	// next, add all the substrates to the microenvironment
	// build the initial conditions and Dirichlet conditions as we go

	// find the first substrate
	pugi::xml_node node1 = node.child("variable"); // xml_find_node( node , "variable" );
	node = node1;
	index_t i = 0;

	while (node)
	{
		// define density
		{
			// get the name and units
			std::string name = node.attribute("name").value();
			std::string units = node.attribute("units").value();

			// get the diffusion and decay parameters
			node1 = xml_find_node(node, "physical_parameter_set");

			auto diffusion_coefficient = xml_get_double_value(node1, "diffusion_coefficient");
			auto decay_rate = xml_get_double_value(node1, "decay_rate");

			// now, get the initial value
			node1 = xml_find_node(node, "initial_condition");
			auto initial_condition = xml_get_my_double_value(node1);

			// add the substrate
			m_builder_.add_density(name, units, diffusion_coefficient, decay_rate, initial_condition);
		}

		// now, get the Dirichlet value
		node1 = xml_find_node(node, "Dirichlet_boundary_condition");
		auto dirichlet_default_value = xml_get_my_double_value(node1);

		// now, decide whether or not to enable it
		auto dirichlet_default_condition = node1.attribute("enabled").as_bool();

		// now figure out finer-grained controls
		point_t<bool, 3> mins_condition = { dirichlet_default_condition, dirichlet_default_condition,
											dirichlet_default_condition };
		point_t<bool, 3> maxs_condition = { dirichlet_default_condition, dirichlet_default_condition,
											dirichlet_default_condition };
		point_t<real_t, 3> mins_value = { dirichlet_default_value, dirichlet_default_value, dirichlet_default_value };
		point_t<real_t, 3> maxs_value = { dirichlet_default_value, dirichlet_default_value, dirichlet_default_value };

		node1 = node.child("Dirichlet_options");
		if (node1)
		{
			// xmin, xmax, ymin, ymax, zmin, zmax, interior
			pugi::xml_node node2 = node1.child("boundary_value");

			while (node2)
			{
				// which boundary?
				std::string boundary_ID = node2.attribute("ID").value();

				index_t dim_idx;

				if (boundary_ID.front() == 'x')
					dim_idx = 0;
				else if (boundary_ID.front() == 'y')
					dim_idx = 1;
				else if (boundary_ID.front() == 'z')
					dim_idx = 2;
				else
				{
					std::cout << "Warning: unknown boundary ID " << boundary_ID << std::endl;
					continue;
				}

				if (boundary_ID.substr(1) == "min")
				{
					mins_condition[dim_idx] = node2.attribute("enabled").as_bool();
					mins_value[dim_idx] = xml_get_my_double_value(node2);
				}
				else if (boundary_ID.substr(1) == "max")
				{
					maxs_condition[dim_idx] = node2.attribute("enabled").as_bool();
					maxs_value[dim_idx] = xml_get_my_double_value(node2);
				}
				else
				{
					std::cout << "Warning: unknown boundary ID " << boundary_ID << std::endl;
					continue;
				}

				node2 = node2.next_sibling("boundary_value");
			}
		}

		m_builder_.add_boundary_dirichlet_conditions(i, mins_value, maxs_value, mins_condition, maxs_condition);

		// move on to the next variable (if any!)
		node = node.next_sibling("variable");
		i++;
	}

	// std::cout << activated_Dirichlet_boundary_detected << std::endl;
	// std::cout << "dc? " << default_microenvironment_options.outer_Dirichlet_conditions << std::endl;

	// now, get the options
	node = xml_find_node(config_root_, "microenvironment_setup");
	node = xml_find_node(node, "options");

	// calculate gradients?
	// default_microenvironment_options.calculate_gradients = xml_get_bool_value(node, "calculate_gradients");

	// track internalized substrates in each agent?
	if (xml_get_bool_value(node, "track_internalized_substrates_in_each_agent"))
		m_builder_.do_compute_internalized_substrates();

	return m_builder_;
}

void builder::load_signals() { setup_signal_behavior_dictionaries(get_environment()); }

void builder::load_rules()
{
	setup_cell_rules(get_settings(), config_root_, get_environment());
	get_environment().rules_enabled = settings_->rules_enabled;
};

void builder::peek_cell_definitions()
{
	if (cell_definition_names_.size())
		return;

	// find the start of cell definitions
	pugi::xml_node node = config_root_.child("cell_definitions");

	// find the first cell definition
	node = node.child("cell_definition");

	cell_definition_names_.push_back("default");

	while (node)
	{
		index_t ID = node.attribute("ID").as_int();
		std::string type_name = node.attribute("name").value();

		std::cout << "Pre-processing type " << ID << " named " << type_name << std::endl;

		if (ID != 0 && type_name != "default")
		{
			cell_definition_names_.push_back(type_name);
		}
		else
		{
			cell_definition_names_[0] = type_name;
		}

		node = node.next_sibling("cell_definition");
	}
}

microenvironment& builder::get_microenvironment()
{
	if (!m_)
		m_.emplace(m_builder_.build());
	return *m_;
}

environment& builder::get_environment()
{
	if (e_)
		return *e_;

	peek_cell_definitions();
	e_.emplace(get_microenvironment(), cell_definition_names_.size());

	pugi::xml_node node_options;

	node_options = xml_find_node(config_root_, "options");
	if (node_options)
	{
		bool settings = xml_get_bool_value(node_options, "virtual_wall_at_domain_edge");
		if (settings)
		{
			std::cout << "virtual_wall_at_domain_edge: enabled" << std::endl;
			e_->virtual_wall_at_domain_edges = true;
		}
	}

	auto overall_node = xml_find_node(config_root_, "overall");

	// check to see if dt is specified in overall options
	// if so, set from XML

	auto search_result = xml_find_node(overall_node, "dt_mechanics");
	if (search_result)
	{
		e_->mechanics_time_step = xml_get_my_double_value(search_result);
	}

	search_result = xml_find_node(overall_node, "dt_phenotype");
	if (search_result)
	{
		e_->phenotype_time_step = xml_get_my_double_value(search_result);
	}

	e_->rules_enabled = get_settings().rules_enabled;
	e_->automated_spring_adhesion = !get_settings().disable_automated_spring_adhesions;

	return *e_;
}

cell_definition& builder::get_default_cell_definition()
{
	auto& e = get_environment();

	if (e.cell_definitions.empty())
		::initialize_default_cell_definition(*e_);

	return e.cell_defaults();
}

std::vector<cell_definition>& builder::get_cell_definitions()
{
	// first, let's pre-build the map.
	// prebuild_cell_definition_index_maps();

	get_default_cell_definition();

	pugi::xml_node node = config_root_.child("cell_definitions");

	node = node.child("cell_definition");

	while (node)
	{
		std::cout << "Processing " << node.attribute("name").value() << " ... " << std::endl;

		construct_single_cell_definition(node);
		// build_cell_definitions_maps();

		node = node.next_sibling("cell_definition");
	}

	return e_->cell_definitions;
}

extern Cycle_Model Ki67_advanced, Ki67_basic, live, apoptosis, necrosis, cycling_quiescent, flow_cytometry_cycle_model,
	flow_cytometry_separated_cycle_model;

void builder::construct_single_cell_definition(const pugi::xml_node& cd_node)
{
	cell_definition* pCD;

	// if this is not "default" then create a new one
	if (std::string(cd_node.attribute("name").value()) != "default"
		&& std::string(cd_node.attribute("ID").value()) != "0")
	{
		e_->cell_definitions.emplace_back(*e_);
		pCD = &e_->cell_definitions.back();
	}
	else
	{
		pCD = &e_->cell_definitions.front();
	}

	// set the name
	pCD->name = cd_node.attribute("name").value();

	// set the ID
	if (cd_node.attribute("ID"))
	{
		pCD->type = cd_node.attribute("ID").as_int();
	}
	else
	{
		pCD->type = -1;
	}

	// get the parent definition (if any)
	cell_definition* pParent = NULL;
	if (cd_node.attribute("parent_type"))
	{
		std::string parent_name = cd_node.attribute("parent_type").value();
		auto it = std::find_if(get_environment().cell_definitions.begin(), get_environment().cell_definitions.end(),
							   [&](const cell_definition& cd) { return cd.name == parent_name; });

		if (it == get_environment().cell_definitions.end())
		{
			throw std::runtime_error("Error: unknown parent type " + parent_name + " in cell definition " + pCD->name);
		}
		pParent = &(*it);
	}

	// if it's not the default and no parent stated, inherit from default
	bool used_default = false;
	if (pParent == nullptr)
	{
		used_default = true;
		pParent = &get_environment().cell_definitions.front();
	}

	// if we found something to inherit from, then do it!
	std::cout << "\tCopying from type " << pParent->name << " ... " << std::endl;
	pCD->inherit_from(*pParent);

	/* bugfix on April 24, 2022 */
	// If we copied from cell_defaults and also wrote
	// more properties to cell defaults, we have messed up
	// some of the rates that are assumed to start at zero.
	// So, let's overwrite with zeros.
	if (used_default)
	{
		pugi::xml_node node_options = xml_find_node(config_root_, "options");
		bool disable_bugfix = false;
		if (node_options)
		{
			xml_get_bool_value(node_options, "legacy_cell_defaults_copy");
		}

		if (disable_bugfix == false)
		{
			// motility
			pCD->phenotype.motility.set_defaults();

			// secretion
			pCD->phenotype.secretion.set_defaults();

			// interaction
			pCD->phenotype.cell_interactions.set_defaults();

			// transformation
			pCD->phenotype.cell_transformations.set_defaults();
		}
		else
		{
			std::cout << "Warning! You have disabled a bugfix on cell definition inheritance" << std::endl
					  << "\tBe VERY careful that you have manually specified every parameter value" << std::endl
					  << "\tfor every cell cell definition. Set legacy_cell_defaults_copy to false" << std::endl
					  << "\tin the options section of your parameter file to re-enable the bug fix. " << std::endl
					  << std::endl
					  << "\tSome good news: if you used the model builder, this has never affected your results."
					  << std::endl;
		}
	}

	// figure out if this ought to be 2D
	if (e_->m.mesh.dims == 2)
	{
		std::cout << "Note: setting cell definition to 2D based on microenvironment domain settings ... " << std::endl;
		pCD->functions.set_orientation = [](cell& cell, real_t) {
			cell.state.orientation()[0] = 0;
			cell.state.orientation()[1] = 0;
		};
		pCD->phenotype.geometry.polarity() = 1.0;
	}

	// make sure phenotype.secretions are correctly sized

	pugi::xml_node node = cd_node.child("phenotype");

	// set up the cell cycle
	// make sure the standard cycle models are defined
	create_standard_cycle_and_death_models();

	node = node.child("cycle");
	if (node)
	{
		index_t model; // = node.attribute("code").as_int() ;
		if (std::string(node.attribute("code").as_string()).size() > 0)
		{
			model = node.attribute("code").as_int();
		}
		else
		{
			model = find_cycle_model_code(node.attribute("name").as_string());
		}
		if (model < 0)
		{
			throw std::runtime_error("Error. Unable to identify cycle model "
									 + std::string(node.attribute("name").value()) + " ("
									 + std::string(node.attribute("code").value()) + ")");
		}

		// set the model
		// switch( model )   // do not use a switch stmt to avoid compile errors related to "static const index_t"
		// on various compilers
		if (model == constants::advanced_Ki67_cycle_model)
		{
			pCD->phenotype.cycle.sync_to_cycle_model(Ki67_advanced);
		}
		else if (model == constants::basic_Ki67_cycle_model)
		{
			pCD->phenotype.cycle.sync_to_cycle_model(Ki67_basic);
		}
		else if (model == constants::flow_cytometry_cycle_model)
		{
			pCD->phenotype.cycle.sync_to_cycle_model(flow_cytometry_cycle_model);
		}
		else if (model == constants::live_apoptotic_cycle_model) // ?
		{
			pCD->phenotype.cycle.sync_to_cycle_model(live);
			std::cout << "Warning: live_apoptotic_cycle_model not directly supported." << std::endl
					  << "         Substituting live cells model. Set death rates=0." << std::endl;
		}
		else if (model == constants::total_cells_cycle_model)
		{
			pCD->phenotype.cycle.sync_to_cycle_model(live);
			std::cout << "Warning: total_cells_cycle_model not directly supported." << std::endl
					  << "         Substituting live cells model. Set death rates=0." << std::endl;
		}
		else if (model == constants::live_cells_cycle_model)
		{
			pCD->phenotype.cycle.sync_to_cycle_model(live);
		}
		else if (model == constants::flow_cytometry_separated_cycle_model)
		{
			pCD->phenotype.cycle.sync_to_cycle_model(flow_cytometry_separated_cycle_model);
		}
		else if (model == constants::cycling_quiescent_model)
		{
			pCD->phenotype.cycle.sync_to_cycle_model(cycling_quiescent);
		}
		else
		{
			throw std::runtime_error("Error. Unknown cycle model " + std::to_string(model));
		}

		// now, if we inherited from another cell, AND
		// if that parent type has the same cylce model,
		// then overwrite with their transition rates

		if (pParent != NULL)
		{
			if (pCD->phenotype.cycle.model().code == pParent->phenotype.cycle.model().code)
			{
				std::cout << "copying data ... " << std::endl;
				std::cout << pParent->name << " to " << pCD->name << std::endl;
				pCD->phenotype.cycle.data = pParent->phenotype.cycle.data;
			}
		}

		// set the transition rates
		if (node.child("phase_transition_rates"))
		{
			node = node.child("phase_transition_rates");
		}
		if (node.child("transition_rates"))
		{
			node = node.child("transition_rates");
			std::cout << "Warning: " << node.name() << " is deprecated. Use cycle.phase_transition_rates." << std::endl;
		}
		if (node)
		{
			node = node.child("rate");
			while (node)
			{
				// which rate
				index_t start = node.attribute("start_index").as_int();
				index_t end = node.attribute("end_index").as_int();
				// fixed duration?
				bool fixed = false;
				if (node.attribute("fixed_duration"))
				{
					fixed = node.attribute("fixed_duration").as_bool();
				}
				// actual value of transition rate
				real_t value = xml_get_my_double_value(node);

				// set the transition rate
				pCD->phenotype.cycle.data.transition_rate(start, end) = value;
				// set it to fixed / non-fixed
				pCD->phenotype.cycle.model().phase_link(start, end).fixed_duration = fixed;

				node = node.next_sibling("rate");
			}
		}

		node = cd_node.child("phenotype");
		node = node.child("cycle");
		// Check for phase durations (as an alternative to transition rates)
		if (node.child("phase_durations"))
		{
			node = node.child("phase_durations");
		}
		if (node)
		{
			node = node.child("duration");
			while (node)
			{
				// which duration?
				index_t start = node.attribute("index").as_int();
				// fixed duration?
				bool fixed = false;
				if (node.attribute("fixed_duration"))
				{
					fixed = node.attribute("fixed_duration").as_bool();
				}
				// actual value of the duration
				real_t value = xml_get_my_double_value(node);

				// set the transition rate
				pCD->phenotype.cycle.data.exit_rate(start) = 1.0 / (value + 1e-16);
				// set it to fixed / non-fixed
				pCD->phenotype.cycle.model().phase_links[start][0].fixed_duration = fixed;

				node = node.next_sibling("duration");
			}
		}
	}

	// here's what it ***should*** do:
	// parse the model, get its code
	// look for that model
	// if the model is not yet there, then add it
	// otherwise, modify properties of that model

	// set up the death models
	//	index_t death_model_index = 0;
	node = cd_node.child("phenotype");
	node = node.child("death");
	if (node)
	{
		pugi::xml_node model_node = node.child("model");
		while (model_node)
		{
			node = model_node;

			index_t model; // = node.attribute("code").as_int() ;
			if (std::string(node.attribute("code").as_string()).size() > 0)
			{
				model = node.attribute("code").as_int();
			}
			else
			{
				model = find_cycle_model_code(node.attribute("name").as_string());
			}
			if (model < 0)
			{
				throw std::runtime_error("Error. Unable to identify death model "
										 + std::string(node.attribute("name").value()) + " ("
										 + std::string(node.attribute("code").value()) + ")");
			}


			// check: is that death model already there?

			Death* pD = &(pCD->phenotype.death);
			index_t death_index = pD->find_death_model_index(model);
			bool death_model_already_exists = false;
			if ((index_t)pD->rates.size() > death_index)
			{
				if (pD->models[death_index]->code == model)
				{
					death_model_already_exists = true;
				}
			}

			// add the death model and its death rate

			if (node.child("death_rate"))
			{
				node = node.child("death_rate");
			}
			if (node.child("rate"))
			{
				node = node.child("rate");
				std::cout << "Warning: " << node.name() << " is deprecated. Use death.model.death_rate." << std::endl;
			}
			real_t rate = xml_get_my_double_value(node);
			node = node.parent();

			// get death model parameters

			Death_Parameters death_params;
			// if there is a parent and we already found this model,
			// start with the inherited parameters
			if (death_model_already_exists && pParent != NULL)
			{
				death_params = pParent->phenotype.death.parameters[death_index];
			}

			if (node.child("parameters"))
			{
				node = node.child("parameters");

				// only read these parameters if they are specified.

				pugi::xml_node node_temp = node.child("unlysed_fluid_change_rate");
				if (node_temp)
				{
					death_params.unlysed_fluid_change_rate = xml_get_my_double_value(node_temp);
				}

				node_temp = node.child("lysed_fluid_change_rate");
				if (node_temp)
				{
					death_params.lysed_fluid_change_rate = xml_get_my_double_value(node_temp);
				}

				node_temp = node.child("cytoplasmic_biomass_change_rate");
				if (node_temp)
				{
					death_params.cytoplasmic_biomass_change_rate = xml_get_my_double_value(node_temp);
				}

				node_temp = node.child("nuclear_biomass_change_rate");
				if (node_temp)
				{
					death_params.nuclear_biomass_change_rate = xml_get_my_double_value(node_temp);
				}

				node_temp = node.child("calcification_rate");
				if (node_temp)
				{
					death_params.calcification_rate = xml_get_my_double_value(node_temp);
				}

				node_temp = node.child("relative_rupture_volume");
				if (node_temp)
				{
					death_params.relative_rupture_volume = xml_get_my_double_value(node_temp);
				}

				node_temp = node.child("lysed_fluid_change_rate");
				if (node_temp)
				{
					death_params.lysed_fluid_change_rate = xml_get_my_double_value(node_temp);
				}

				//			death_params.time_units =
				//				get_string_attribute_value( node, "unlysed_fluid_change_rate", "units" );

				node = node.parent();
			}

			// set the model
			// if the model already exists, just overwrite the parameters
			if (model == constants::apoptosis_death_model)
			{
				//					pCD->phenotype.death.add_death_model( rate , &apoptosis , apoptosis_parameters
				//);
				if (death_model_already_exists == false)
				{
					pCD->phenotype.death.add_death_model(rate, &apoptosis, death_params);
					death_index = pD->find_death_model_index(model);
				}
				else
				{
					pCD->phenotype.death.parameters[death_index] = death_params;
					pCD->phenotype.death.rates[death_index] = rate;
				}
			}
			else if (model == constants::necrosis_death_model)
			{
				// set necrosis parameters
				//					pCD->phenotype.death.add_death_model( rate , &necrosis , necrosis_parameters );
				if (death_model_already_exists == false)
				{
					pCD->phenotype.death.add_death_model(rate, &necrosis, death_params);
					death_index = pD->find_death_model_index(model);
				}
				else
				{
					pCD->phenotype.death.parameters[death_index] = death_params;
					pCD->phenotype.death.rates[death_index] = rate;
				}
			}
			else if (model == constants::autophagy_death_model)
			{
				std::cout << "Warning: autophagy_death_model not yet supported." << std::endl
						  << "         Skipping this model." << std::endl;
			}
			else
			{
				throw std::invalid_argument("Unknown death model");
			}

			// now get transition rates within the death model
			// set the rates
			// node = node.child( "transition_rates" );


			pugi::xml_node node_death_transitions = node.child("phase_transition_rates");
			if (node.child("transition_rates"))
			{
				node_death_transitions = node.child("transition_rates");
				std::cout << "Warning: " << node_death_transitions.name()
						  << " is deprecated. Use death.model.phase_transition_rates." << std::endl;
			}


			if (node_death_transitions)
			{
				pugi::xml_node node1 = node_death_transitions.child("rate");
				while (node1)
				{
					// which rate
					index_t start = node1.attribute("start_index").as_int();
					index_t end = node1.attribute("end_index").as_int();
					// fixed duration?
					bool fixed = false;
					if (node1.attribute("fixed_duration"))
					{
						fixed = node1.attribute("fixed_duration").as_bool();
					}
					// actual value of transition rate
					real_t value = xml_get_my_double_value(node1);

					// set the transition rate
					pCD->phenotype.death.models[death_index]->transition_rate(start, end) = value;
					// set it to fixed / non-fixed
					pCD->phenotype.death.models[death_index]->phase_link(start, end).fixed_duration = fixed;

					node1 = node1.next_sibling("rate");
				}
			}

			if (node.child("phase_durations"))
			{
				node = node.child("phase_durations"); // phase durations
				node = node.child("duration");		  // duration
				while (node)
				{
					// which duration?
					index_t start = node.attribute("index").as_int();
					// fixed duration?
					bool fixed = false;
					if (node.attribute("fixed_duration"))
					{
						fixed = node.attribute("fixed_duration").as_bool();
					}
					// actual value of the duration
					real_t value = xml_get_my_double_value(node);

					// set the transition rate
					pCD->phenotype.death.models[death_index]->data.exit_rate(start) = 1.0 / (value + 1e-16);
					// set it to fixed / non-fixed
					pCD->phenotype.death.models[death_index]->phase_links[start][0].fixed_duration = fixed;

					node = node.next_sibling("duration");
				}


				/*
						if( node.child( "phase_durations" ) )
						{ node = node.child( "phase_durations" ); }
						if( node )
						{
							node = node.child( "duration");
							while( node )
							{
								// which duration?
								index_t start = node.attribute("index").as_int();
								// fixed duration?
								bool fixed = false;
								if( node.attribute( "fixed_duration" ) )
								{ fixed = node.attribute("fixed_duration").as_bool(); }
								// actual value of the duration
								real_t value = xml_get_my_double_value( node );

								// set the transition rate
								pCD->phenotype.cycle.data.exit_rate(start) = 1.0 / (value+1e-16);
								// set it to fixed / non-fixed
								pCD->phenotype.cycle.model().phase_links[start][0].fixed_duration = fixed;

								node = node.next_sibling( "duration" );
							}

						}

				*/


				node = node.parent(); // phase_durations
				node = node.parent(); // model
			}

			// node = node.parent();

			model_node = model_node.next_sibling("model");
			//			death_model_index++;
		}
	}

	// volume
	node = cd_node.child("phenotype");
	node = node.child("volume");
	if (node)
	{
		volume_t* pV = &(pCD->phenotype.volume);

		pugi::xml_node node_vol = node.child("total");
		if (node_vol)
		{
			pV->total() = xml_get_my_double_value(node_vol);
		}

		node_vol = node.child("fluid_fraction");
		if (node_vol)
		{
			pV->fluid_fraction() = xml_get_my_double_value(node_vol);
		}

		node_vol = node.child("nuclear");
		if (node_vol)
		{
			pV->nuclear() = xml_get_my_double_value(node_vol);
		}

		node_vol = node.child("fluid_change_rate");
		if (node_vol)
		{
			pV->fluid_change_rate() = xml_get_my_double_value(node_vol);
		}

		node_vol = node.child("cytoplasmic_biomass_change_rate");
		if (node_vol)
		{
			pV->cytoplasmic_biomass_change_rate() = xml_get_my_double_value(node_vol);
		}

		node_vol = node.child("nuclear_biomass_change_rate");
		if (node_vol)
		{
			pV->nuclear_biomass_change_rate() = xml_get_my_double_value(node_vol);
		}

		node_vol = node.child("calcified_fraction");
		if (node_vol)
		{
			pV->calcified_fraction() = xml_get_my_double_value(node_vol);
		}

		node_vol = node.child("calcification_rate");
		if (node_vol)
		{
			pV->calcification_rate() = xml_get_my_double_value(node_vol);
		}

		node_vol = node.child("relative_rupture_volume");
		if (node_vol)
		{
			pV->relative_rupture_volume() = xml_get_my_double_value(node_vol);
		}

		// set all the parameters to be self-consistent

		pV->fluid() = pV->fluid_fraction() * pV->total();
		pV->solid() = pV->total() - pV->fluid();

		pV->nuclear_fluid() = pV->fluid_fraction() * pV->nuclear();
		pV->nuclear_solid() = pV->nuclear() - pV->nuclear_fluid();

		pV->cytoplasmic() = pV->total() - pV->nuclear();
		pV->cytoplasmic_fluid() = pV->fluid_fraction() * pV->cytoplasmic();
		pV->cytoplasmic_solid() = pV->cytoplasmic() - pV->cytoplasmic_fluid();


		pV->target_solid_cytoplasmic() = pV->cytoplasmic_solid();
		pV->target_solid_nuclear() = pV->nuclear_solid();
		pV->target_fluid_fraction() = pV->fluid_fraction();

		pV->cytoplasmic_to_nuclear_ratio() = pV->cytoplasmic() / (1e-16 + pV->nuclear());
		pV->target_cytoplasmic_to_nuclear_ratio() = pV->cytoplasmic_to_nuclear_ratio();

		pV->rupture_volume() = pV->relative_rupture_volume() * pV->total(); // in volume units

		// update the geometry (radius, etc.) for consistency

		pCD->phenotype.geometry.update();
	}

	// mechanics
	node = cd_node.child("phenotype");
	node = node.child("mechanics");
	if (node)
	{
		mechanics_t* pM = &(pCD->phenotype.mechanics);

		pugi::xml_node node_mech = node.child("cell_cell_adhesion_strength");
		if (node_mech)
		{
			pM->cell_cell_adhesion_strength() = xml_get_my_double_value(node_mech);
		}

		node_mech = node.child("cell_BM_adhesion_strength");
		if (node_mech)
		{
			pM->cell_BM_adhesion_strength() = xml_get_my_double_value(node_mech);
		}

		node_mech = node.child("cell_cell_repulsion_strength");
		if (node_mech)
		{
			pM->cell_cell_repulsion_strength() = xml_get_my_double_value(node_mech);
		}

		node_mech = node.child("cell_BM_repulsion_strength");
		if (node_mech)
		{
			pM->cell_BM_repulsion_strength() = xml_get_my_double_value(node_mech);
		}

		node_mech = node.child("relative_maximum_adhesion_distance");
		if (node_mech)
		{
			pM->relative_maximum_adhesion_distance() = xml_get_my_double_value(node_mech);
		}

		// cell adhesion affinities
		node_mech = node.child("cell_adhesion_affinities");
		if (node_mech)
		{
			node_mech = node_mech.child("cell_adhesion_affinity");
			while (node_mech)
			{
				std::string target = node_mech.attribute("name").value();
				real_t value = xml_get_my_double_value(node_mech);

				// find the target
				// if found, assign taht affinity
				index_t ind = find_cell_definition_index(target);
				if (ind > -1)
				{
					pM->cell_adhesion_affinities()[ind] = value;
				}
				else
				{
					throw std::invalid_argument("cell adhesion affinity target not found:" + target);
				}

				node_mech = node_mech.next_sibling("cell_adhesion_affinity");
			}
		}

		node_mech = node.child("options");
		if (node_mech)
		{
			pugi::xml_node node_mech1 = node_mech.child("set_relative_equilibrium_distance");
			if (node_mech1)
			{
				if (node_mech1.attribute("enabled").as_bool())
				{
					real_t temp = xml_get_my_double_value(node_mech1);
					pM->set_relative_equilibrium_distance(temp);
				}
			}

			node_mech1 = node_mech.child("set_absolute_equilibrium_distance");
			if (node_mech1)
			{
				if (node_mech1.attribute("enabled").as_bool())
				{
					real_t temp = xml_get_my_double_value(node_mech1);
					pM->set_absolute_equilibrium_distance(pCD->phenotype.geometry.radius(), temp);
				}
			}
		}

		node_mech = node.child("attachment_elastic_constant");
		if (node_mech)
		{
			pM->attachment_elastic_constant() = xml_get_my_double_value(node_mech);
		}
		std::cout << "  --------- attachment_elastic_constant = " << pM->attachment_elastic_constant() << std::endl;

		node_mech = node.child("attachment_rate");
		if (node_mech)
		{
			pM->attachment_rate() = xml_get_my_double_value(node_mech);
		}

		node_mech = node.child("detachment_rate");
		if (node_mech)
		{
			pM->detachment_rate() = xml_get_my_double_value(node_mech);
		}
	}

	// motility
	node = cd_node.child("phenotype");
	node = node.child("motility");
	if (node)
	{
		motility_t* pMot = &(pCD->phenotype.motility);
		decltype(chemotaxis_function)* motility_f = nullptr;

		pugi::xml_node node_mot = node.child("speed");
		if (node_mot)
		{
			pMot->migration_speed() = xml_get_my_double_value(node_mot);
		}

		node_mot = node.child("migration_bias");
		if (node_mot)
		{
			pMot->migration_bias() = xml_get_my_double_value(node_mot);
		}

		node_mot = node.child("persistence_time");
		if (node_mot)
		{
			pMot->persistence_time() = xml_get_my_double_value(node_mot);
		}

		node_mot = node.child("options");
		if (node_mot)
		{
			// enable motility?
			pugi::xml_node node_mot1 = node_mot.child("enabled");
			if (node_mot1)
			{
				pMot->is_motile() = xml_get_my_bool_value(node_mot1);
			}

			node_mot1 = node_mot.child("use_2D");

			if ((e_->m.mesh.dims == 2 && xml_get_my_bool_value(node_mot1) == false)
				|| (e_->m.mesh.dims == 3 && xml_get_my_bool_value(node_mot1) == true))
			{
				throw std::invalid_argument(
					"Error: use_2D must be set to true for 2D simulations and false for 3D simulations.");
			}

			// automated chemotaxis setup
			node_mot1 = node_mot.child("chemotaxis");
			if (node_mot1)
			{
				// enabled? if so, set the standard chemotaxis function
				if (xml_get_bool_value(node_mot1, "enabled"))
				{
					motility_f = chemotaxis_function;
					pMot->update_migration_bias_direction() = chemotaxis_function;
				}

				// search for the right chemo index

				std::string substrate_name = xml_get_string_value(node_mot1, "substrate");
				pMot->chemotaxis_index() = e_->m.find_substrate_index(substrate_name);
				if (pMot->chemotaxis_index() < 0)
				{
					throw std::invalid_argument(
						"Error: parsing phenotype:motility:options:chemotaxis:  invalid substrate " + substrate_name);
				}

				// set the direction

				pMot->chemotaxis_direction() = xml_get_int_value(node_mot1, "direction");

				// std::cout << pMot->chemotaxis_direction << " * grad( " << actual_name << " )" << std::endl;
			}

			// automated advanced chemotaxis setup
			node_mot1 = node_mot.child("advanced_chemotaxis");
			if (node_mot1)
			{
				// enabled? if so, set the standard chemotaxis function
				if (xml_get_bool_value(node_mot1, "enabled"))
				{
					if (motility_f == chemotaxis_function)
					{
						std::cout << "Warning: when processing motility for " << pCD->name << " cells: " << std::endl
								  << "\tBoth chemotaxis and advanced_chemotaxis are enabled." << std::endl
								  << "\tThe settings for advanced_chemotaxis override those of chemotaxis."
								  << std::endl;
					}
					pMot->update_migration_bias_direction() = advanced_chemotaxis_function<false>;
					motility_f = advanced_chemotaxis_function<false>;
					if (xml_get_bool_value(node_mot1, "normalize_each_gradient"))
					{
						motility_f = advanced_chemotaxis_function<true>;
						pMot->update_migration_bias_direction() = advanced_chemotaxis_function<true>;
					}
				}

				// now process the chemotactic sensitivities

				pugi::xml_node node_cs = node_mot1.child("chemotactic_sensitivities");
				if (node_cs)
				{
					node_cs = node_cs.child("chemotactic_sensitivity");

					while (node_cs)
					{
						std::string substrate_name = node_cs.attribute("substrate").value();
						index_t index = e_->m.find_substrate_index(substrate_name);
						std::string actual_name = "";
						if (index == -1)
						{
							throw std::invalid_argument("Error: parsing "
														"phenotype:motility:options:advanced_chemotaxis:chemotactic_"
														"sensitivities:  invalid substrate "
														+ substrate_name);
						}
						pCD->phenotype.motility.chemotactic_sensitivities()[index] = xml_get_my_double_value(node_cs);
						node_cs = node_cs.next_sibling("chemotactic_sensitivity");
					}
				}
				else
				{
					std::cout << "Warning: when processing motility for " << pCD->name << " cells: " << std::endl
							  << "\tAdvanced chemotaxis requries chemotactic_sensitivities." << std::endl
							  << "\tBut you have none. Your migration bias will be the zero vector." << std::endl;
				}
			}
		}

		// display summary for diagnostic help
		if (motility_f == chemotaxis_function && pMot->is_motile() == true)
		{
			std::cout << "Cells of type " << pCD->name << " use standard chemotaxis: " << std::endl
					  << "\t d_bias (before normalization) = " << pMot->chemotaxis_direction() << " * grad("
					  << e_->m.substrates_names[pMot->chemotaxis_index()] << ")" << std::endl;
		}

		if (motility_f == advanced_chemotaxis_function<false> && pMot->is_motile() == true)
		{
			index_t number_of_substrates = e_->m.substrates_count;

			std::cout << "Cells of type " << pCD->name << " use advanced chemotaxis: " << std::endl
					  << "\t d_bias (before normalization) = " << pMot->chemotactic_sensitivities()[0] << " * grad("
					  << e_->m.substrates_names[0] << ")";

			for (index_t n = 1; n < number_of_substrates; n++)
			{
				std::cout << " + " << pMot->chemotactic_sensitivities()[n] << " * grad(" << e_->m.substrates_names[n]
						  << ")";
			}
			std::cout << std::endl;
		}

		if (motility_f == advanced_chemotaxis_function<true> && pMot->is_motile() == true)
		{
			index_t number_of_substrates = e_->m.substrates_count;

			std::cout << "Cells of type " << pCD->name << " use normalized advanced chemotaxis: " << std::endl
					  << "\t d_bias (before normalization) = " << pMot->chemotactic_sensitivities()[0] << " * grad("
					  << e_->m.substrates_names[0] << ")"
					  << " / ||grad(" << e_->m.substrates_names[0] << ")||";

			for (index_t n = 1; n < number_of_substrates; n++)
			{
				std::cout << " + " << pMot->chemotactic_sensitivities()[n] << " * grad(" << e_->m.substrates_names[n]
						  << ")"
						  << " / ||grad(" << e_->m.substrates_names[n] << ")||";
			}
			std::cout << std::endl;
		}
	}

	// secretion

	node = cd_node.child("phenotype");
	node = node.child("secretion");
	if (node)
	{
		secretion_t* pS = &(pCD->phenotype.secretion);

		// find the first substrate
		pugi::xml_node node_sec = node.child("substrate");
		while (node_sec)
		{
			// which substrate?

			std::string substrate_name = node_sec.attribute("name").value();
			index_t index = e_->m.find_substrate_index(substrate_name);

			// error check
			if (index == -1)
			{
				throw std::invalid_argument("Error: parsing phenotype:secretion:substrate:  invalid substrate "
											+ substrate_name);
			}

			// secretion rate
			pugi::xml_node node_sec1 = node_sec.child("secretion_rate");
			if (node_sec1)
			{
				pS->secretion_rates()[index] = xml_get_my_double_value(node_sec1);
			}

			// secretion target
			node_sec1 = node_sec.child("secretion_target");
			if (node_sec1)
			{
				pS->saturation_densities()[index] = xml_get_my_double_value(node_sec1);
			}

			// uptake rate
			node_sec1 = node_sec.child("uptake_rate");
			if (node_sec1)
			{
				pS->uptake_rates()[index] = xml_get_my_double_value(node_sec1);
			}

			// net export rate
			node_sec1 = node_sec.child("net_export_rate");
			if (node_sec1)
			{
				pS->net_export_rates()[index] = xml_get_my_double_value(node_sec1);
			}

			node_sec = node_sec.next_sibling("substrate");
		}
	}

	// cell interactions

	node = cd_node.child("phenotype");
	node = node.child("cell_interactions");
	if (node)
	{
		interactions_t* pCI = &(pCD->phenotype.cell_interactions);

		// dead_phagocytosis_rate
		pugi::xml_node node_dpr = node.child("dead_phagocytosis_rate");
		pCI->dead_phagocytosis_rate() = xml_get_my_double_value(node_dpr);

		// live phagocytosis rates
		pugi::xml_node node_lpcr = node.child("live_phagocytosis_rates");
		if (node_lpcr)
		{
			node_lpcr = node_lpcr.child("phagocytosis_rate");
		}
		while (node_lpcr)
		{
			// get the name of the target cell type
			std::string target_name = node_lpcr.attribute("name").value();
			// now find its index
			index_t index = find_cell_definition_index(target_name);
			// safety first!
			if (index == -1)
			{
				throw std::invalid_argument(
					"Error: parsing phenotype:cell_interactions:live_phagocytosis_rates:  invalid cell type "
					+ target_name);
			}

			// if the target is found, set the appropriate rate
			pCI->live_phagocytosis_rates()[index] = xml_get_my_double_value(node_lpcr);

			node_lpcr = node_lpcr.next_sibling("phagocytosis_rate");
		}

		// effector attack rates
		pugi::xml_node node_ar = node.child("attack_rates");
		if (node_ar)
		{
			node_ar = node_ar.child("attack_rate");
		}
		while (node_ar)
		{
			// get the name of the target cell type
			std::string target_name = node_ar.attribute("name").value();
			// now find its index
			index_t index = find_cell_definition_index(target_name);
			// safety first!
			if (index == -1)
			{
				throw std::invalid_argument(
					"Error: parsing phenotype:cell_interactions:attack_rates:  invalid cell type " + target_name);
			}
			pCI->attack_rates()[index] = xml_get_my_double_value(node_ar);

			node_ar = node_ar.next_sibling("attack_rate");
		}

		// damage_rate
		pugi::xml_node node_dr = node.child("damage_rate");
		pCI->damage_rate() = xml_get_my_double_value(node_dr);

		// fusion_rates
		pugi::xml_node node_fr = node.child("fusion_rates");
		if (node_fr)
		{
			node_fr = node_fr.child("fusion_rate");
		}
		while (node_fr)
		{
			// get the name of the target cell type
			std::string target_name = node_fr.attribute("name").value();
			// now find its index
			index_t index = find_cell_definition_index(target_name);
			// safety first!
			if (index == -1)
			{
				throw std::invalid_argument(
					"Error: parsing phenotype:cell_interactions:fusion_rate:  invalid cell type " + target_name);
			}
			pCI->fusion_rates()[index] = xml_get_my_double_value(node_ar);

			node_fr = node_fr.next_sibling("fusion_rate");
		}
	}

	// cell_transformations>
	//            <transformation_rate

	node = cd_node.child("phenotype");
	node = node.child("cell_transformations");
	if (node)
	{
		transformations_t* pCT = &(pCD->phenotype.cell_transformations);

		// transformation rates
		pugi::xml_node node_tr = node.child("transformation_rates");
		if (node_tr)
		{
			node_tr = node_tr.child("transformation_rate");
		}
		while (node_tr)
		{
			// get the name of the target cell type
			std::string target_name = node_tr.attribute("name").value();
			// now find its index
			index_t index = find_cell_definition_index(target_name);
			// safety first!
			if (index == -1)
			{
				throw std::invalid_argument(
					"Error: parsing phenotype:cell_interactions:transformation_rate:  invalid cell type "
					+ target_name);
			}
			// if the target is found, set the appropriate rate

			real_t transformation_rate = xml_get_my_double_value(node_tr);
			if (target_name == pCD->name && transformation_rate > 1e-16)
			{
				std::cout << "Warning: When processing the " << pCD->name << " cell definition: " << std::endl
						  << "\tTransformation from " << pCD->name << " to " << target_name << " is not allowed."
						  << std::endl
						  << "\tIgnoring this cell transformation rate!" << std::endl
						  << std::endl;
			}
			else
			{
				pCT->transformation_rates()[index] = transformation_rate;
			}

			node_tr = node_tr.next_sibling("transformation_rate");
		}
	}

	// intracellular
	node = cd_node.child("phenotype");
	node = node.child("intracellular");
	if (node)
	{
		std::string model_type = node.attribute("type").value();


#ifdef ADDON_PHYSIBOSS
		if (model_type == "maboss")
		{
			// If it has already be copied
			if (pParent != NULL && pParent->phenotype.intracellular != NULL)
			{
				pCD->phenotype.intracellular->initialize_intracellular_from_pugixml(node);

				// Otherwise we need to create a new one
			}
			else
			{
				MaBoSSIntracellular* pIntra = new MaBoSSIntracellular(node);
				pCD->phenotype.intracellular = pIntra->getIntracellularModel();
			}
		}
#endif

#ifdef ADDON_ROADRUNNER
		if (model_type == "roadrunner")
		{
			// If it has already be copied
			if (pParent != NULL && pParent->phenotype.intracellular != NULL)
			{
				// std::cout << "------ " << __FUNCTION__ << ": copying another\n";
				pCD->phenotype.intracellular->initialize_intracellular_from_pugixml(node);
			}
			// Otherwise we need to create a new one
			else
			{
				std::cout << "\n------ " << __FUNCTION__ << ": creating new RoadRunnerIntracellular\n";
				RoadRunnerIntracellular* pIntra = new RoadRunnerIntracellular(node);
				pCD->phenotype.intracellular = pIntra->getIntracellularModel();
				pCD->phenotype.intracellular->validate_PhysiCell_tokens(pCD->phenotype);
				pCD->phenotype.intracellular->validate_SBML_species();
			}
		}
#endif

#ifdef ADDON_PHYSIDFBA
		if (model_type == "dfba")
		{
			// If it has already be copied
			if (pParent != NULL && pParent->phenotype.intracellular != NULL)
			{
				pCD->phenotype.intracellular->initialize_intracellular_from_pugixml(node);
				// Otherwise we need to create a new one
			}
			else
			{
				dFBAIntracellular* pIntra = new dFBAIntracellular(node);
				pCD->phenotype.intracellular = pIntra->getIntracellularModel();
			}
		}
#endif
	}

	// set up custom data
	node = cd_node.child("custom_data");
	pugi::xml_node node1 = node.first_child();
	while (node1)
	{
		// name of teh custom data
		std::string name = xml_get_my_name(node1);

		// units
		std::string units = node1.attribute("units").value();

		// conserved quantity
		bool conserved = node1.attribute("conserved").as_bool();

		// get value(s)
		std::string str_values = xml_get_my_string_value(node1);
		std::vector<real_t> values = csv_to_vector(str_values.c_str());

		// add variable if cell defaults
		// if the custom data is not yet found, add it
		// first, try scalar
		if (values.size() == 1)
		{
			// find the variable
			index_t n = pCD->custom_data.find_variable_index(name);
			// if it exists, overwrite
			if (n > -1)
			{
				pCD->custom_data.variables[n].value = values[0];
				pCD->custom_data.variables[n].conserved_quantity = conserved;
			}
			// otherwise, add
			else
			{
				auto index = pCD->custom_data.add_variable(name, units, values[0]);
				pCD->custom_data.variables[index].conserved_quantity = conserved;
			}

			n = pCD->custom_data.find_variable_index(name);
		}
		// or vector
		else
		{
			// find the variable
			index_t n = pCD->custom_data.find_vector_variable_index(name);
			// if it exists, overwrite
			if (n > -1)
			{
				pCD->custom_data.vector_variables[n].value = values;
				pCD->custom_data.vector_variables[n].conserved_quantity = conserved;
			}
			// otherwise, add
			else
			{
				auto index = pCD->custom_data.add_vector_variable(name, units, values);
				pCD->custom_data.vector_variables[index].conserved_quantity = conserved;
			}

			n = pCD->custom_data.find_vector_variable_index(name);
		}

		// set conserved attribute

		node1 = node1.next_sibling();
	}
}

environment builder::build_environment()
{
	get_cell_definitions();

	load_rules();

	load_signals();

	random::instance().set_seed(get_parameters().ints("random_seed"));

	e_->display_info();

	return std::move(e_.value());
};
