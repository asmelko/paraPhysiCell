#pragma once

#include "common_solver.h"

namespace physicell {

class position_solver : public common_solver
{
public:
	static void update_cell_forces(environment& e);

	static void update_cell_forces_new(environment& e);

	static void update_cell_neighbors(environment& e);

	static void update_motility(environment& e);

	static void update_basement_membrane_interactions(environment& e);

	static void update_spring_attachments(environment& e);

	static void update_positions(environment& e);
	
	static void update_positions_new(environment& e);
};

} // namespace physicell
