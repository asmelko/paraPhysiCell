#pragma once

#include "common_solver.h"

namespace physicell {

class containers_solver : public common_solver
{
public:
	static void update_mechanics_mesh(environment& e);

	static void update_cell_container_for_mechanics(environment& e);
	static void update_cell_container_for_phenotype(environment& e, biofvm::cell_solver& s);
};

void remove_springs(biofvm::index_t to_remove, std::vector<biofvm::index_t>* __restrict__ springs);

} // namespace physicell
