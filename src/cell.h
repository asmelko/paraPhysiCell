#pragma once

#include <cstdint>
#include <vector>

#include <BioFVM/agent.h>

#include "cell_data.h"
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

	const environment& e() const;

	biofvm::real_t* velocity();
	biofvm::index_t& cell_definition_index();
	std::uint8_t& is_movable();
	cell_state_flag& flag();

	void set_default(cell_definition& def, biofvm::index_t def_index);
	void convert(cell_definition& def, biofvm::index_t def_index);
	void copy_from(cell& source);

	void divide(cell& new_cell);

	void assign_orientation();

	void attach_cell(cell& cell);
	void detach_cell(cell& cell);

	static void attach_cells(cell& lhs, cell& rhs);
	static void detach_cells(cell& lhs, cell& rhs);

	void flag_for_removal();
	void flag_for_division();
};

} // namespace physicell
