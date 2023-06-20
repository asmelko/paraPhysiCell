#pragma once

#include <cstdint>
#include <vector>

#include <BioFVM/agent.h>

#include "phenotype.h"

namespace physicell {

struct cell_data;
class cell_container;
class position_solver;

class cell : public biofvm::agent
{
	friend cell_container;
	friend position_solver;

	cell_data& data_;

	std::vector<cell*> neighbors_;
	std::vector<cell*> spring_attached_cells_;

public:
	cell(biofvm::agent_id_t id, cell_data& data, biofvm::index_t index);

	phenotype_t phenotype;

	biofvm::real_t* velocity();
	biofvm::index_t& cell_definition_index();
	biofvm::real_t& simple_pressure();
	std::uint8_t& is_movable();

	std::vector<cell*>& neighbors();
	std::vector<cell*>& spring_attached_cells();

	biofvm::index_t index() const { return index_; }
};

} // namespace physicell
