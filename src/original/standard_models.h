#pragma once

#include "../environment.h"

namespace physicell {

void standard_cell_transformations(cell& pCell, biofvm::real_t time_step);

void advance_bundled_phenotype_functions(environment& e);

void chemotaxis_function(cell& cell);

template <bool do_normalize>
void advanced_chemotaxis_function(cell& cell);

void standard_volume_update_function(cell& cell);

void standard_elastic_contract_function(cell& lhs, cell& rhs);

void evaluate_interactions(environment& e);

void update_cell_and_death_parameters_O2_based(cell& cell, biofvm::real_t dt);

} // namespace physicell
