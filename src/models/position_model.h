#pragma once

namespace physicell {

class environment;

class position_model
{
public:
	virtual void update_cell_neighbors(environment& e) = 0;

	virtual void update_cell_positions(environment& e) = 0;
};

} // namespace physicell
