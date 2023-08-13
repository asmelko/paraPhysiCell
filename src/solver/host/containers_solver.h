#pragma once

#include <mutex>
#include <shared_mutex>

#include "common_solver.h"

namespace physicell {

class containers_solver : public common_solver
{
	std::mutex removal_mtx_;
	std::shared_mutex division_mtx_;

public:
	static void update_mechanics_mesh(environment& e);

	void update_cell_container_for_mechanics(environment& e);
	void update_cell_container_for_phenotype(environment& e, biofvm::cell_solver& s);
};

void remove_attached(biofvm::index_t to_remove, std::vector<biofvm::index_t>* __restrict__ springs);

} // namespace physicell
