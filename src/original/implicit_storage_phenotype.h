#pragma once

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

#include "../cell_functions.h"
#include "../data_storage_phenotype.h"

namespace physicell {

class Phase
{
public:
	biofvm::index_t index; // an internal index for the cycle model
	biofvm::index_t code;  // a global identifier code
	std::string name;

	bool division_at_phase_exit; // does this phase trigger division?
	bool removal_at_phase_exit;	 // does this phase trigger removal?

	cell_func_t<void> entry_function;

	Phase(); // done
};

class Phase_Link
{
public:
	biofvm::index_t start_phase_index;
	biofvm::index_t end_phase_index;

	bool fixed_duration;

	cell_func_t<bool> arrest_function;
	// return true if arrested, false if not

	cell_func_t<void> exit_function;
	// function to be excecuted when completing the phase transition

	Phase_Link(); // done
};

class Cycle_Model;

class Cycle_Data
{
private:
	// this maps the end_phase_index to the link index in each
	// phase_links[i]
	// So, index_inverse_map[i][j] = k, corresponds to
	// phases[i], phase_links[i][k] (which links from phase i to phase j)
	// transition_rates[i][k] (the transition rate from phase i to phase j)
	std::vector<std::unordered_map<biofvm::index_t, biofvm::index_t>> inverse_index_maps;

public:
	Cycle_Model* pCycle_Model;

	std::string time_units;

	std::vector<std::vector<biofvm::real_t>> transition_rates;

	biofvm::index_t current_phase_index;
	biofvm::real_t elapsed_time_in_phase;

	Cycle_Data(); // done

	// return current phase (by reference)
	Phase& current_phase(); // done

	// make the data structures consistent with the corresponding cell cycle model
	void sync_to_cycle_model(); // done

	// access the transition rate from phase i to phase j (by reference)
	biofvm::real_t& transition_rate(biofvm::index_t start_phase_index, biofvm::index_t end_phase_index); // done

	biofvm::real_t& exit_rate(biofvm::index_t phase_index); // This returns the first transition rate out of
															// phase # phase_index. It is only relevant if the phase has
															// only one phase link (true for many cycle models).
};

class Cycle_Model
{
private:
	// this maps the end_phase_index to the link index in each
	// phase_links[i]
	// So, index_inverse_map[i][j] = k, corresponds to
	// phases[i], phase_links[i][k] (which links from phase i to phase j)
	// transition_rates[i][k] (the transition rate from phase i to phase j)
	std::vector<std::unordered_map<biofvm::index_t, biofvm::index_t>> inverse_index_maps;

public:
	std::string name;
	biofvm::index_t code;

	std::vector<Phase> phases;
	std::vector<std::vector<Phase_Link>> phase_links;

	biofvm::index_t default_phase_index;

	Cycle_Data data; // this will be copied to individual cell agents

	Cycle_Model();

	void advance_model(cell& Cell, biofvm::real_t dt); // done

	biofvm::index_t add_phase(biofvm::index_t code, std::string name); // done

	biofvm::index_t add_phase_link(biofvm::index_t start_index, biofvm::index_t end_index,
								   cell_func_t<bool> arrest_function); // done
	biofvm::index_t add_phase_link(biofvm::index_t start_index, biofvm::index_t end_index, biofvm::real_t rate,
								   cell_func_t<bool> arrest_function); // done

	biofvm::index_t find_phase_index(biofvm::index_t code); // done
	biofvm::index_t find_phase_index(std::string name);		// done

	biofvm::real_t& transition_rate(biofvm::index_t start_index, biofvm::index_t end_index); // done
	Phase_Link& phase_link(biofvm::index_t start_index, biofvm::index_t end_index);			 // done

	std::ostream& display(std::ostream& os); // done
};

class Cycle
{
public:
	Cycle_Model* pCycle_Model;
	Cycle_Data data;

	Cycle(); // done

	void advance_cycle(cell& pCell, biofvm::real_t dt); // done

	Cycle_Model& model();					// done
	Phase& current_phase();					// done
	biofvm::index_t& current_phase_index(); // done

	void sync_to_cycle_model(Cycle_Model& cm); // done

	void copy(Cycle& dest);
};

class Death_Parameters
{
public:
	std::string time_units;

	biofvm::real_t unlysed_fluid_change_rate;
	biofvm::real_t lysed_fluid_change_rate;

	biofvm::real_t cytoplasmic_biomass_change_rate;
	biofvm::real_t nuclear_biomass_change_rate;

	biofvm::real_t calcification_rate;

	biofvm::real_t relative_rupture_volume;

	Death_Parameters(); // done
};

class Death : public phenotype_data_storage
{
public:
	std::vector<biofvm::real_t> rates;
	std::vector<Cycle_Model*> models;
	std::vector<Death_Parameters> parameters;

	std::uint8_t& dead();
	biofvm::index_t current_death_model_index;

	Death(cell_data& data, const biofvm::index_t& index); // done

	biofvm::index_t add_death_model(biofvm::real_t rate, Cycle_Model* pModel); // done
	biofvm::index_t add_death_model(biofvm::real_t rate, Cycle_Model* pModel,
									Death_Parameters& death_parameters); // done

	biofvm::index_t find_death_model_index(biofvm::index_t code); // done
	biofvm::index_t find_death_model_index(std::string name);	  // done

	bool check_for_death(biofvm::real_t dt);			   // done
	void trigger_death(biofvm::index_t death_model_index); // done

	Cycle_Model& current_model();			// done
	Death_Parameters& current_parameters(); // done

	// ease of access
	biofvm::real_t& apoptosis_rate();
	biofvm::real_t& necrosis_rate();

	void copy(Death& dest);
};

} // namespace physicell
