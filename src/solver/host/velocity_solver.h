#pragma once

#include "../../cell_data.h"

namespace physicell {

class velocity_solver
{
public:
	static void solve(cell_data& data);
};

} // namespace physicell
