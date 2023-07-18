#pragma once

#include "common_solver.h"

namespace physicell {

class mechanics_solver : public common_solver
{
public:
	static void update_mechanics_mesh(environment& e);

	static void update_cell_container(environment& e);
};

} // namespace physicell
