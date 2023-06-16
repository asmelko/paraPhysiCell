#pragma once

#include "../../environment.h"

namespace physicell {

class mechanics_solver
{
public:
	static void update_mechanics_mesh(environment& e);

	static void update_cell_neighbors(environment& e);
};

} // namespace physicell
