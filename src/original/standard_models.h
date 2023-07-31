#pragma once

#include "../environment.h"

namespace physicell {

void standard_cell_transformations(cell& pCell, biofvm::real_t time_step);

void advance_bundled_phenotype_functions(environment& e);

void chemotaxis_function(cell& cell);

template <bool do_normalize>
void advanced_chemotaxis_function(cell& cell);

} // namespace physicell
