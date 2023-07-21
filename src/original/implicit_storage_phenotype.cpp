#include "implicit_storage_phenotype.h"

#include <cstddef>
#include <iostream>

#include "../cell.h"
#include "../cell_data.h"
#include "../random.h"
#include "constants.h"

using namespace biofvm;
using namespace physicell;

Phase::Phase()
{
	index = 0;
	code = 0;
	name = "unnamed";

	division_at_phase_exit = false;
	removal_at_phase_exit = false;

	entry_function = NULL;
	return;
}

Phase_Link::Phase_Link()
{
	start_phase_index = 0;
	end_phase_index = 0;

	fixed_duration = false;

	arrest_function = NULL;
	exit_function = NULL;
	return;
}

Cycle_Data::Cycle_Data()
{
	inverse_index_maps.resize(0);

	pCycle_Model = NULL;

	time_units = "min";

	transition_rates.resize(0);

	current_phase_index = 0;
	elapsed_time_in_phase = 0.0;
	return;
}

void Cycle_Data::sync_to_cycle_model(void)
{
	// make sure the inverse map is the right size
	biofvm::index_t n = pCycle_Model->phases.size();
	inverse_index_maps.resize(n);

	// sync the inverse map to the cell cycle model by
	// querying the phase_links

	transition_rates.resize(n);

	// also make sure the transition_rates[] are the right size

	for (std::size_t i = 0; i < pCycle_Model->phase_links.size(); i++)
	{
		inverse_index_maps[i].clear();
		for (std::size_t j = 0; j < pCycle_Model->phase_links[i].size(); j++)
		{
			inverse_index_maps[i][pCycle_Model->phase_links[i][j].end_phase_index] = j;
			transition_rates[i].resize(pCycle_Model->phase_links[i].size());
		}
	}

	return;
}

biofvm::real_t& Cycle_Data::transition_rate(biofvm::index_t start_phase_index, biofvm::index_t end_phase_index)
{
	return transition_rates[start_phase_index][inverse_index_maps[start_phase_index][end_phase_index]];
}

biofvm::real_t& Cycle_Data::exit_rate(biofvm::index_t phase_index) { return transition_rates[phase_index][0]; }

Cycle_Model::Cycle_Model()
{
	inverse_index_maps.resize(0);

	name = "unnamed";

	phases.resize(0);
	phase_links.resize(0);

	data.pCycle_Model = this;

	code = constants::custom_cycle_model;
	default_phase_index = 0;

	return;
}

biofvm::index_t Cycle_Model::add_phase(biofvm::index_t code, std::string name)
{
	biofvm::index_t n = phases.size();

	// resize the data structures
	phases.resize(n + 1);
	phase_links.resize(n + 1);
	phase_links[n].resize(0);

	inverse_index_maps.resize(n + 1);
	inverse_index_maps[n].clear();

	// update phase n
	phases[n].code = code;
	phases[n].index = n;

	phases[n].name.assign(name);

	// make sure the cycle_data is also correctly sized

	data.sync_to_cycle_model();

	return n;
}

biofvm::index_t Cycle_Model::add_phase_link(biofvm::index_t start_index, biofvm::index_t end_index,
											cell_func_t<bool>& arrest_function)
{
	// first, resize the phase links
	biofvm::index_t n = phase_links[start_index].size();
	phase_links[start_index].resize(n + 1);

	// now, update the new phase links
	phase_links[start_index][n].start_phase_index = start_index;
	phase_links[start_index][n].end_phase_index = end_index;
	phase_links[start_index][n].arrest_function = arrest_function;

	// now, update the inverse index map
	inverse_index_maps[start_index][end_index] = n;

	// lastly, make sure the transition rates are the right size;

	data.sync_to_cycle_model();

	return n;
}

biofvm::index_t Cycle_Model::add_phase_link(biofvm::index_t start_index, biofvm::index_t end_index, biofvm::real_t rate,
											cell_func_t<bool>& arrest_function)
{
	biofvm::index_t n = add_phase_link(start_index, end_index, arrest_function);
	data.transition_rate(start_index, end_index) = rate;
	return n;
}

biofvm::index_t Cycle_Model::find_phase_index(biofvm::index_t code)
{
	for (std::size_t i = 0; i < phases.size(); i++)
	{
		if (phases[i].code == code)
		{
			return i;
		}
	}
	return 0;
}

biofvm::index_t Cycle_Model::find_phase_index(std::string name)
{
	for (std::size_t i = 0; i < phases.size(); i++)
	{
		if (phases[i].name == name)
		{
			return i;
		}
	}
	return 0;
}

std::ostream& Cycle_Model::display(std::ostream& os)
{
	os << "Cycle Model: " << name << " (PhysiCell code: " << code << ")" << std::endl;
	os << "Phases and links: (* denotes phase with cell division)" << std::endl;
	for (std::size_t i = 0; i < phases.size(); i++)
	{
		os << "Phase " << i << " (" << phases[i].name << ") ";

		if (phases[i].division_at_phase_exit)
		{
			os << "*";
		}
		os << " links to: " << std::endl;
		for (std::size_t k = 0; k < phase_links[i].size(); k++)
		{
			std::size_t j = phase_links[i][k].end_phase_index;
			os << "\tPhase " << j << " (" << phases[j].name << ") with rate " << data.transition_rate(i, j) << " "
			   << data.time_units << "^-1; " << std::endl;
		}
		os << std::endl;
	}

	return os;
}

biofvm::real_t& Cycle_Model::transition_rate(biofvm::index_t start_index, biofvm::index_t end_index)
{
	return data.transition_rate(start_index, end_index);
}

Phase_Link& Cycle_Model::phase_link(biofvm::index_t start_index, biofvm::index_t end_index)
{
	return phase_links[start_index][inverse_index_maps[start_index][end_index]];
}

void Cycle_Model::advance_model(cell& cell, biofvm::real_t dt)
{
	biofvm::index_t i = cell.phenotype.cycle.data.current_phase_index;

	cell.phenotype.cycle.data.elapsed_time_in_phase += dt;

	// Evaluate each linked phase:
	// advance to that phase IF probabiltiy is in the range,
	// and if the arrest function (if any) is false

	biofvm::index_t j;
	for (std::size_t k = 0; k < phase_links[i].size(); k++)
	{
		j = phase_links[i][k].end_phase_index;

		// check for arrest. If arrested, skip to the next transition
		bool transition_arrested = false;
		if (phase_links[i][k].arrest_function)
		{
			transition_arrested = phase_links[i][k].arrest_function(cell, dt);
		}
		if (!transition_arrested)
		{
			// check to see if we should transition
			bool continue_transition = false;
			if (phase_links[i][k].fixed_duration)
			{
				if (cell.phenotype.cycle.data.elapsed_time_in_phase
					> 1.0 / cell.phenotype.cycle.data.transition_rates[i][k])
				{
					continue_transition = true;
				}
			}
			else
			{
				biofvm::real_t prob = cell.phenotype.cycle.data.transition_rates[i][k] * dt;
				if (random::instance().uniform() < prob)
				{
					continue_transition = true;
				}
			}

			// if we should transition, check if we're not supposed to divide or die

			if (continue_transition)
			{
				// if the phase transition has an exit function, execute it
				if (phase_links[i][k].exit_function)
				{
					phase_links[i][k].exit_function(cell, dt);
				}

				// check if division or removal are required
				if (phases[i].division_at_phase_exit)
				{
					cell.divide();
				}
				if (phases[i].removal_at_phase_exit)
				{
					cell.remove();
					return;
				}
				// move to the next phase, and reset the elapsed time
				cell.phenotype.cycle.data.current_phase_index = j;
				cell.phenotype.cycle.data.elapsed_time_in_phase = 0.0;

				// if the new phase has an entry function, execute it
				if (phases[j].entry_function)
				{
					phases[j].entry_function(cell, dt);
				}

				return;
			}
		}
	}

	return;
}

Phase& Cycle_Data::current_phase(void) { return pCycle_Model->phases[current_phase_index]; }

Death_Parameters::Death_Parameters()
{
	time_units = "min";

	// reference values: MCF-7 (1/min)
	unlysed_fluid_change_rate = 3.0 / 60.0; // apoptosis
	lysed_fluid_change_rate = 0.05 / 60.0;	// lysed necrotic cell

	cytoplasmic_biomass_change_rate = 1.0 / 60.0; // apoptosis
	nuclear_biomass_change_rate = 0.35 / 60.0;	  // apoptosis

	calcification_rate = 0.0; // 0.0042 for necrotic cells

	relative_rupture_volume = 2.0;

	return;
}

Death::Death(cell_data& data, index_t index) : phenotype_data_storage(data, index)
{
	rates.resize(0);
	models.resize(0);
	parameters.resize(0);

	dead() = false;
	current_death_model_index = 0;

	return;
}

std::uint8_t& Death::dead() { return data_.deaths.dead[index_]; }

biofvm::index_t Death::add_death_model(biofvm::real_t rate, Cycle_Model* pModel)
{
	rates.push_back(rate);
	models.push_back(pModel);

	parameters.resize(rates.size());

	return rates.size() - 1;
}

biofvm::index_t Death::add_death_model(biofvm::real_t rate, Cycle_Model* pModel, Death_Parameters& death_parameters)
{
	rates.push_back(rate);
	models.push_back(pModel);
	parameters.push_back(death_parameters);

	return rates.size() - 1;
}

biofvm::index_t Death::find_death_model_index(biofvm::index_t code)
{
	for (std::size_t i = 0; i < models.size(); i++)
	{
		if (models[i]->code == code)
		{
			return i;
		}
	}
	return 0;
}

biofvm::index_t Death::find_death_model_index(std::string name)
{
	for (std::size_t i = 0; i < models.size(); i++)
	{
		if (models[i]->name == name)
		{
			return i;
		}
	}
	return 0;
}

bool Death::check_for_death(biofvm::real_t dt)
{
	// If the cell is already dead, exit.
	if (dead() == true)
	{
		return false;
	}

	// If the cell is alive, evaluate all the
	// death rates for each registered death type.
	std::size_t i = 0;
	while (!dead() && i < rates.size())
	{
		if (random::instance().uniform() < rates[i] * dt)
		{
			// update the Death data structure
			dead() = true;
			current_death_model_index = i;

			// and set the cycle model to this death model

			return dead();
		}
		i++;
	}

	return dead();
}

void Death::trigger_death(biofvm::index_t death_model_index)
{
	dead() = true;
	current_death_model_index = death_model_index;

	/*
		// if so, change the cycle model to the current death model
		phenotype.cycle.sync_to_cycle_model( phenotype.death.current_model() );

		// also, turn off motility.

		phenotype.motility.is_motile = false;
		phenotype.motility.motility_vector.assign( 3, 0.0 );
		functions.update_migration_bias = NULL;

		// turn off secretion, and reduce uptake by a factor of 10
		phenotype.secretion.set_all_secretion_to_zero();
		phenotype.secretion.scale_all_uptake_by_factor( 0.10 );

		// make sure to run the death entry function
		if( phenotype.cycle.current_phase().entry_function )
		{
			phenotype.cycle.current_phase().entry_function( this, phenotype, dt_ );
		}
	*/

	return;
}

Cycle_Model& Death::current_model(void) { return *models[current_death_model_index]; }

biofvm::real_t& Death::apoptosis_rate(void)
{
	static biofvm::index_t nApoptosis = find_death_model_index(constants::apoptosis_death_model);
	return rates[nApoptosis];
}

biofvm::real_t& Death::necrosis_rate(void)
{
	static biofvm::index_t nNecrosis = find_death_model_index(constants::necrosis_death_model);
	return rates[nNecrosis];
}


Cycle::Cycle()
{
	pCycle_Model = NULL;
	return;
}

void Cycle::advance_cycle(cell& cell, biofvm::real_t dt)
{
	pCycle_Model->advance_model(cell, dt);
	return;
}

Cycle_Model& Cycle::model(void) { return *pCycle_Model; }

Phase& Cycle::current_phase(void) { return data.current_phase(); }

biofvm::index_t& Cycle::current_phase_index(void) { return data.current_phase_index; }

void Cycle::sync_to_cycle_model(Cycle_Model& cm)
{
	pCycle_Model = &cm;
	data = cm.data;
	return;
}

Death_Parameters& Death::current_parameters(void) { return parameters[current_death_model_index]; }

void Cycle::copy(Cycle& dest) { dest = *this; }

void Death::copy(Death& dest)
{
	dest.rates = rates;
	dest.models = models;
	dest.parameters = parameters;
	dest.dead() = dead();
	dest.current_death_model_index = current_death_model_index;
}
