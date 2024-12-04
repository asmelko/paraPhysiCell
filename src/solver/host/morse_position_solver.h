#pragma once

#include "common_solver.h"

namespace physicell {

class morse_position_solver : public common_solver
{
public:
	static void update_cell_forces(environment& e);

	static void update_motility(environment& e);

	static void update_positions(environment& e);
};

} // namespace physicell
