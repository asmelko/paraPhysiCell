#pragma once

#include <vector>

#include <BioFVM/agent.h>

#include "phenotype.h"

namespace physicell {

struct cell_data;
class cell_container;
class velocity_solver;

class cell : public biofvm::agent
{
	friend cell_container;
	friend velocity_solver;

	cell_data& data_;

	std::vector<cell*> neighbors_;

public:
	cell(biofvm::agent_id_t id, cell_data& data, biofvm::index_t index);

	phenotype_t phenotype;

	biofvm::real_t* velocity();
	biofvm::index_t& cell_definition_index();
	biofvm::real_t& simple_pressure();

	std::vector<cell*>& neighbors();
};

} // namespace physicell
