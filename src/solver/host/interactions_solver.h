#pragma once

#include "common_solver.h"

namespace physicell {

class interactions_solver : public common_solver
{
public:
	static void update_cell_cell_interactions(environment& e);
};

void update_geometry(biofvm::index_t i, biofvm::real_t* __restrict__ radius,
					 biofvm::real_t* __restrict__ nuclear_radius, biofvm::real_t* __restrict__ surface_area,
					 const biofvm::real_t* __restrict__ total_volume,
					 const biofvm::real_t* __restrict__ nuclear_volume);

} // namespace physicell
