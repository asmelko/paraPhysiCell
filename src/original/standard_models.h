#pragma once

#include "../environment.h"

namespace physicell {

// standard entry function for the cycle models

void standard_Ki67_positive_phase_entry_function(cell& cell, biofvm::real_t dt); // done
void standard_Ki67_negative_phase_entry_function(cell& cell, biofvm::real_t dt); // done
void standard_live_phase_entry_function(cell& cell, biofvm::real_t dt);			 // done

void G1_phase_entry_function(cell& cell, biofvm::real_t dt);
void G0_phase_entry_function(cell& cell, biofvm::real_t dt);
void S_phase_entry_function(cell& cell, biofvm::real_t dt); // done

void standard_apoptosis_entry_function(cell& cell, biofvm::real_t dt); // done
void standard_necrosis_entry_function(cell& cell, biofvm::real_t dt);  // done
void standard_lysis_entry_function(cell& cell, biofvm::real_t dt);	   // done

bool standard_necrosis_arrest_function(cell& cell, biofvm::real_t dt); // done

// create standard models

bool create_standard_cell_cycle_models();	   // done
bool create_standard_cell_death_models();	   // done
bool create_standard_cycle_and_death_models(); // done

void initialize_default_cell_definition(environment& e);

// standard volume functions

void standard_volume_update_function(cell& cell);

// THESE ARE OPTIMIZED:
// standard mechanics functions
// bounary avoidance functions

// standard o2-based phenotype changes

void update_cell_and_death_parameters_O2_based(cell& cell);

// standard motility functions

void chemotaxis_function(cell& cell);

template <bool do_normalize>
void advanced_chemotaxis_function(cell& cell);

// automated cell phagocytosis, attack, and fusion

void evaluate_interactions(environment& e);

// standard contact functions

void standard_elastic_contract_function(cell& lhs, cell& rhs);

// standard phenotype update functions

void advance_bundled_phenotype_functions(environment& e);
void standard_cell_transformations(cell& pCell, biofvm::real_t time_step);

} // namespace physicell
