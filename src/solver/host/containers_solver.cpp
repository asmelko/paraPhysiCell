#include "containers_solver.h"

#include <algorithm>

#include "BioFVM/solver.h"
#include "solver_helper.h"

using namespace biofvm;
using namespace physicell;

void containers_solver::initialize(environment& e)
{
	cells_in_voxels_sizes_ = std::make_unique<std::atomic<index_t>[]>(e.mechanics_mesh.voxel_count());
}

void containers_solver::update_mechanics_mesh(environment& e)
{
#pragma omp for
	for (std::size_t i = 0; i < e.mechanics_mesh.voxel_count(); i++)
	{
		e.cells_in_mechanics_voxels[i].clear();
		cells_in_voxels_sizes_[i].store(0, std::memory_order_relaxed);
	}

	const auto& data = get_cell_data(e);

	// first we count how many cells are in each voxel
#pragma omp for
	for (index_t i = 0; i < data.agents_count; i++)
	{
		auto mech_pos =
			get_mesh_position(data.agent_data.positions.data() + i * e.mechanics_mesh.dims, e.mechanics_mesh);

		cells_in_voxels_sizes_[get_mesh_index(mech_pos, e.mechanics_mesh)].fetch_add(1, std::memory_order_relaxed);
	}

	// second we allocate memory for each voxel
#pragma omp for
	for (std::size_t i = 0; i < e.mechanics_mesh.voxel_count(); i++)
	{
		e.cells_in_mechanics_voxels[i].resize(cells_in_voxels_sizes_[i].load(std::memory_order_relaxed));
	}

	// third we assign cells to voxels
#pragma omp for
	for (index_t i = 0; i < data.agents_count; i++)
	{
		auto mech_pos =
			get_mesh_position(data.agent_data.positions.data() + i * e.mechanics_mesh.dims, e.mechanics_mesh);

		auto mech_idx = get_mesh_index(mech_pos, e.mechanics_mesh);

		auto in_voxel_index = cells_in_voxels_sizes_[mech_idx].fetch_sub(1, std::memory_order_relaxed) - 1;

		e.cells_in_mechanics_voxels[mech_idx][in_voxel_index] = i;
	}
}

void physicell::remove_attached(index_t to_remove, std::vector<index_t>* __restrict__ attached)
{
	for (auto spring : attached[to_remove])
	{
		auto it = std::find_if(attached[spring].begin(), attached[spring].end(),
							   [to_remove](index_t i) { return i == to_remove; });

		*it = attached[spring].back();
		attached[spring].pop_back();
	}
}

void rename_attached(index_t old_index, index_t new_index, std::vector<index_t>* __restrict__ attached)
{
	for (auto spring : attached[old_index])
	{
		auto it = std::find_if(attached[spring].begin(), attached[spring].end(),
							   [old_index](index_t i) { return i == old_index; });

		*it = new_index;
	}
}

void remove_single(index_t i, const cell_state_flag* __restrict__ flag, const real_t* __restrict__ positions,
				   std::vector<index_t>* __restrict__ springs, std::vector<index_t>* __restrict__ attachced_cells,
				   const cartesian_mesh& mesh, cell_container& container, index_t& counter, std::mutex& m)
{
	while (true)
	{
		bool out_of_bounds = false;

		for (index_t d = 0; d < mesh.dims; d++)
		{
			if (positions[i * mesh.dims + d] < mesh.bounding_box_mins[d]
				|| positions[i * mesh.dims + d] > mesh.bounding_box_maxs[d])
			{
				out_of_bounds = true;
			}
		}

		if (out_of_bounds == false && flag[i] != cell_state_flag::to_remove)
			return;

		{
			std::lock_guard<std::mutex> l(m);

			// this can happen, it is safe because we are not deallocating memory for removed elements
			if (i >= container.data().agents_count)
				return;

			if (flag[i] == cell_state_flag::to_remove)
				counter++;

			remove_attached(i, springs);
			rename_attached(container.data().agents_count - 1, i, springs);

			remove_attached(i, attachced_cells);
			rename_attached(container.data().agents_count - 1, i, attachced_cells);

			container.remove_at(i);

			if (i == container.data().agents_count)
				return;
		}
	}
}

void update_cell_container_internal(const cell_state_flag* __restrict__ flag, const real_t* __restrict__ positions,
									std::vector<index_t>* __restrict__ springs,
									std::vector<index_t>* __restrict__ attached_cells, const cartesian_mesh& mesh,
									cell_container& container, index_t& counter, std::mutex& m)
{
	auto n = container.data().agents_count;
#pragma omp barrier

#pragma omp for
	for (index_t i = 0; i < n; i++)
	{
		remove_single(i, flag, positions, springs, attached_cells, mesh, container, counter, m);
	}
}

void containers_solver::update_cell_container_for_mechanics(environment& e)
{
	auto& data = get_cell_data(e);

	update_cell_container_internal(data.flags.data(), data.agent_data.positions.data(), data.states.springs.data(),
								   data.states.attached_cells.data(), e.m.mesh, e.get_container(), e.deaths_count,
								   removal_mtx_);
}

void containers_solver::update_cell_container_for_phenotype(environment& e, cell_solver& s)
{
	auto& data = get_cell_data(e);
	auto& c = e.get_container();

	const auto n = c.agents().size();
#pragma omp barrier

#pragma omp for
	for (std::size_t i = 0; i < n; i++)
	{
		if (c.agents()[i]->flag() == cell_state_flag::to_divide)
		{
			c.agents()[i]->flag() = cell_state_flag::none;

			cell* cell;
			{
				std::unique_lock<std::shared_mutex> l(division_mtx_);

				cell = c.create();
				e.divisions_count++;
			}

			{
				std::shared_lock<std::shared_mutex> l(division_mtx_);

				c.agents()[i]->divide(*cell);
			}
		}
	}

#pragma omp for
	for (std::size_t i = 0; i < n; i++)
	{
		s.release_internalized_substrates(e.m, i);

		remove_single(i, data.flags.data(), data.agent_data.positions.data(), data.states.springs.data(),
					  data.states.attached_cells.data(), e.m.mesh, c, e.deaths_count, removal_mtx_);
	}
}
