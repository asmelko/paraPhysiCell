#pragma once

#include "cell_data.h"
#include "cell_functions.h"
#include "original/cell_parameters.h"
#include "original/custom_cell_data.h"
#include "phenotype.h"

namespace physicell {

struct environment;

struct cell_definition
{
private:
	cell_data& data_;

public:
	biofvm::index_t index;

	biofvm::index_t type;
	std::string name;

	bool is_movable;

	environment& e;

	Cell_Parameters parameters;
	Custom_Cell_Data custom_data;
	cell_functions functions;
	phenotype_t phenotype;

	cell_definition(environment& e, biofvm::index_t index);

	// TODO: fix so it does not have to be always unique ptr
	cell_definition(cell_definition&&) = delete;

	void inherit_from(cell_definition& def);
};

} // namespace physicell
