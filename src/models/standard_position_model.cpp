#include "standard_position_model.h"

#include "../solver/host/standard_position_solver.h"

using namespace physicell;

void standard_position_model::update_cell_positions(environment& e)
{
	update_cell_forces(e);
	update_motility(e);
	update_basement_membrane_interactions(e);
	update_spring_attachments(e);
	update_positions(e);
}
