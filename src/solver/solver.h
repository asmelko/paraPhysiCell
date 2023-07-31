#pragma once

#define impl host

#if impl == host
#	include "host/position_solver.h"
#	include "host/containers_solver.h"
#	include "host/interactions_solver.h"
#endif

namespace physicell {

class solver
{
public:
	position_solver position;
	containers_solver containers;
	interactions_solver interactions;

	void initialize(environment& e);
};

} // namespace physicell