#pragma once

#include "standard_position_model.h"

namespace physicell {

class morse_position_model : public standard_position_model
{
protected:
	static void update_cell_forces(environment& e);

	static void update_motility(environment& e);

	static void update_positions(environment& e);

public:
	virtual void update_cell_positions(environment& e) override;
};

} // namespace physicell
