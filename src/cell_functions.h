#pragma once

#include <functional>

#include <BioFVM/types.h>

namespace physicell {

class cell;

template <typename T>
using cell_func_t = std::function<T(cell& cell, biofvm::real_t dt)>;

struct cell_functions
{
	cell_func_t<void> volume_update_function = nullptr;

	cell_func_t<void> custom_cell_rule = nullptr;
	cell_func_t<void> update_phenotype = nullptr;

	cell_func_t<void> pre_update_intracellular = nullptr;
	cell_func_t<void> post_update_intracellular = nullptr;

	cell_func_t<void> set_orientation = nullptr;
};

} // namespace physicell
