#pragma once

#include "environment.h"

namespace physicell {

struct Cell_Definition
{
private:
	cell_data data_;

public:
	int type;
	std::string name;

	bool is_movable;

	environment& e;

	Cell_Parameters parameters;
	Custom_Cell_Data custom_data;
	// Cell_Functions functions;
	phenotype_t phenotype;

	Cell_Definition(environment& e);
};

} // namespace physicell
