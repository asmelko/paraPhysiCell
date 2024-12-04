#pragma once

#include "position_model.h"

namespace physicell {

class standard_position_model : public position_model
{
protected:
	static void update_cell_forces(environment& e);

	static void update_motility(environment& e);

	static void update_basement_membrane_interactions(environment& e);

	static void update_spring_attachments(environment& e);

	static void update_positions(environment& e);

public:
	virtual void update_cell_neighbors(environment& e) override;

	virtual void update_cell_positions(environment& e) override;
};

} // namespace physicell
