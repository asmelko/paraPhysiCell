#pragma once

#include "common_solver.h"

namespace physicell {

class interactions_solver : public common_solver
{
public:
	static void update_cell_cell_interactions(environment& e);
};

} // namespace physicell
