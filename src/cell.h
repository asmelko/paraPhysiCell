#pragma once

#include <cstdint>
#include <vector>

#include <BioFVM/agent.h>

#include "cell_definition.h"

namespace physicell {

struct cell_data;

struct cell_state_t : public phenotype_data_storage
{
	cell_state_t(cell_data& data, biofvm::index_t index);

	std::vector<biofvm::index_t>& neighbors();

	std::vector<biofvm::index_t>& spring_attachments();
	std::vector<biofvm::index_t>& attached_cells();

	biofvm::real_t* orientation();
	biofvm::real_t& simple_pressure();
	biofvm::index_t& number_of_nuclei();

	biofvm::real_t& damage();
	biofvm::real_t& total_attack_time();
};

class cell : public biofvm::agent
{
	cell_data& data_;

public:
	cell(biofvm::agent_id_t id, cell_data& data, biofvm::index_t index);

	biofvm::index_t type;
	std::string type_name;

	Cell_Parameters parameters;

	Custom_Cell_Data custom_data;

	cell_functions functions;

	phenotype_t phenotype;

	cell_state_t state;

	biofvm::real_t* velocity();
	biofvm::index_t& cell_definition_index();
	std::uint8_t& is_movable();

	void set_default(cell_definition& def);
	void convert(cell_definition& def);

	void assign_orientation();

	// TODO:
	void remove();
	void divide();
};

} // namespace physicell
