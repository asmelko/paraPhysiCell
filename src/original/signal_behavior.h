#pragma once

#include <string>
#include <vector>

#include <BioFVM/types.h>

namespace physicell {

struct cell_definition;
struct environment;
class cell;

// scales for the signals
extern std::vector<biofvm::real_t> signal_scales;
// easy access to get or set scales
biofvm::real_t& signal_scale(std::string signal_name);		// done
biofvm::real_t& signal_scale(biofvm::index_t signal_index); // done

// create the signal and behavior dictionaries
void setup_signal_behavior_dictionaries(environment& e); // done

// display dictionaries
void display_signal_dictionary(void);	// done
void display_behavior_dictionary(void); // done

void display_signal_dictionary(std::ostream& os);	// done
void display_behavior_dictionary(std::ostream& os); // done

void display_signal_dictionary_with_synonyms(void);	  // done
void display_behavior_dictionary_with_synonyms(void); // done

/* signal functions */

// find index for named signal (returns -1 if not found)
biofvm::index_t find_signal_index(std::string signal_name); // done

// coming soon:
std::vector<biofvm::index_t> find_signal_indices(std::vector<std::string> signal_names); // done

// get the name of a signal index
std::string signal_name(biofvm::index_t i); // done

// create a full signal vector
std::vector<biofvm::real_t> get_signals(cell* pCell, environment& e); // done

// create a signal vector of only the cell contacts
std::vector<biofvm::real_t> get_cell_contact_signals(cell* pCell, environment& e); // done

// create a subset of the signal vector with the supplied indicies
std::vector<biofvm::real_t> get_selected_signals(cell* pCell, std::vector<biofvm::index_t> indices,
												 environment& e);											   // done
std::vector<biofvm::real_t> get_selected_signals(cell* pCell, std::vector<std::string> names, environment& e); // done

// grab a single signal by its index or name
biofvm::real_t get_single_signal(cell* pCell, biofvm::index_t index, environment& e); // done
biofvm::real_t get_single_signal(cell* pCell, std::string name);					  // done

/* behavior functions */

// find index for named behavior / response / parameter (returns -1 if not found)
biofvm::index_t find_parameter_index(std::string response_name); // done
biofvm::index_t find_behavior_index(std::string response_name);	 // done

std::vector<biofvm::index_t> find_behavior_indices(std::vector<std::string> behavior_names); // done

// get the name of a behavior index
std::string behavior_name(biofvm::index_t i); // done

// make a properly sized behavior vector
std::vector<biofvm::real_t> create_empty_behavior_vector(); // done

// write a full behavior vector (phenotype parameters) to the cell
void set_behaviors(cell* pCell, std::vector<biofvm::real_t> parameters); // done

// write a selected set of behavior parameters to the cell
void set_selected_behaviors(cell* pCell, std::vector<biofvm::index_t> indices,
							std::vector<biofvm::real_t> parameters); // done
void set_selected_behaviors(cell* pCell, std::vector<std::string> names,
							std::vector<biofvm::real_t> parameters); // done

// write a single behavior parameter
void set_single_behavior(cell* pCell, biofvm::index_t index, biofvm::real_t parameter); // done
void set_single_behavior(cell* pCell, std::string name, biofvm::real_t parameter);		// done

/* get current behaviors */

// get all current behavior
std::vector<biofvm::real_t> get_behaviors(cell* pCell); // done

// get selected current behavior
std::vector<biofvm::real_t> get_behaviors(cell* pCell, std::vector<biofvm::index_t> indices); // doen
std::vector<biofvm::real_t> get_behaviors(cell* pCell, std::vector<std::string> names);		  // done

// get single current behavior
biofvm::real_t get_single_behavior(cell* pCell, biofvm::index_t index, environment& e); // done
biofvm::real_t get_single_behavior(cell* pCell, std::string name, environment& e);		// done

/* get base behaviors (from cell definition) */

// get all base behaviors (from cell's definition)
std::vector<biofvm::real_t> get_base_behaviors(cell* pCell, environment& e); // done

// get selected base behaviors (from cell's definition)
std::vector<biofvm::real_t> get_base_behaviors(cell* pCell, std::vector<biofvm::index_t> indices,
											   environment& e);												 // done
std::vector<biofvm::real_t> get_base_behaviors(cell* pCell, std::vector<std::string> names, environment& e); // done

// get single base behavior (from cell's definition)
biofvm::real_t get_single_base_behavior(cell* pCell, biofvm::index_t index, environment& e); // done
biofvm::real_t get_single_base_behavior(cell* pCell, std::string name, environment& e);		 // done

biofvm::real_t get_single_base_behavior(cell_definition* pCD, std::string name, environment& e);

}; // namespace physicell
