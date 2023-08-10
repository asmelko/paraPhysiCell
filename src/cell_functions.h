#pragma once

#include <functional>

#include <BioFVM/types.h>

namespace physicell {

class cell;

template <typename T>
using cell_func_t = std::function<T(cell& cell, biofvm::real_t dt)>;

struct cell_functions
{
	std::function<void(cell& cell)> volume_update_function = nullptr;

	std::function<void(cell& cell)> custom_cell_rule = nullptr;
	std::function<void(cell& cell)> update_phenotype = nullptr;

	cell_func_t<void> pre_update_intracellular = nullptr;
	cell_func_t<void> post_update_intracellular = nullptr;

	cell_func_t<void> set_orientation = nullptr;

	std::function<void(cell& lhs, cell& rhs)> contact_function = nullptr;
};

} // namespace physicell
