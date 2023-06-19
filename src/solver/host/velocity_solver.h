#pragma once

#include "common_solver.h"

namespace physicell {

class velocity_solver : public common_solver
{
public:
	static void solve(environment& e);
};

} // namespace physicell
