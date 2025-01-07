#pragma once

#include "common_solver.h"

namespace physicell {

class morse_kv_membrane_solver : public common_solver
{
public:
	static void update_cell_forces(environment& e);

	static void update_cell_neighbors(environment& e);
};

} // namespace physicell
