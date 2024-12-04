#pragma once

#include "host/containers_solver.h"
#include "host/interactions_solver.h"

namespace physicell {

class solver
{
public:
	containers_solver containers;
	interactions_solver interactions;

	void initialize(environment& e);
};

} // namespace physicell
