#include "velocity_solver.h"

#include <cmath>

#include <BioFVM/microenvironment.h>

using namespace biofvm;
using namespace physicell;

constexpr real_t simple_pressure_coefficient = 36.64504274775163; // 1 / (12 * (1 - sqrt(pi/(2*sqrt(3))))^2)

void clear_simple_pressure(real_t* __restrict__ simple_pressure, index_t count)
{
	for (index_t i = 0; i < count; i++)
	{
		simple_pressure[i] = 0;
	}
}

template <index_t dims>
struct velocity_helper
{};

template <>
struct velocity_helper<1>
{
	static constexpr real_t difference_and_distance(const real_t* __restrict__ lhs, const real_t* __restrict__ rhs,
													real_t* __restrict__ difference)
	{
		difference[0] = lhs[0] - rhs[0];

		return std::abs(difference[0]);
	}

	static constexpr void update_velocities(real_t* __restrict__ lhs, real_t* __restrict__ rhs,
											const real_t* __restrict__ difference, const real_t force)
	{
		lhs[0] += force * difference[0];
		rhs[0] -= force * difference[0];
	}
};

template <>
struct velocity_helper<2>
{
	static constexpr real_t difference_and_distance(const real_t* __restrict__ lhs, const real_t* __restrict__ rhs,
													real_t* __restrict__ difference)
	{
		difference[0] = lhs[0] - rhs[0];
		difference[1] = lhs[1] - rhs[1];

		return std::sqrt(difference[0] * difference[0] + difference[1] * difference[1]);
	}

	static constexpr void update_velocities(real_t* __restrict__ lhs, real_t* __restrict__ rhs,
											const real_t* __restrict__ difference, const real_t force)
	{
		lhs[0] += force * difference[0];
		lhs[1] += force * difference[1];

		rhs[0] -= force * difference[0];
		rhs[1] -= force * difference[1];
	}
};

template <>
struct velocity_helper<3>
{
	static constexpr real_t difference_and_distance(const real_t* __restrict__ lhs, const real_t* __restrict__ rhs,
													real_t* __restrict__ difference)
	{
		difference[0] = lhs[0] - rhs[0];
		difference[1] = lhs[1] - rhs[1];
		difference[2] = lhs[2] - rhs[2];

		return std::sqrt(difference[0] * difference[0] + difference[1] * difference[1] + difference[2] * difference[2]);
	}

	static constexpr void update_velocities(real_t* __restrict__ lhs, real_t* __restrict__ rhs,
											const real_t* __restrict__ difference, const real_t force)
	{
		lhs[0] += force * difference[0];
		lhs[1] += force * difference[1];
		lhs[2] += force * difference[2];

		rhs[0] -= force * difference[0];
		rhs[1] -= force * difference[1];
		rhs[2] -= force * difference[2];
	}
};

template <index_t dims>
void solve_pair(index_t lhs, index_t rhs, real_t* __restrict__ velocity, real_t* __restrict__ simple_pressure,
				const real_t* __restrict__ position, const real_t* __restrict__ radius,
				const real_t* __restrict__ cell_cell_repulsion_strength,
				const real_t* __restrict__ cell_cell_adhesion_strength,
				const real_t* __restrict__ relative_maximum_adhesion_distance, real_t lhs_cell_adhesion_affinity,
				real_t rhs_cell_adhesion_affinity, std::unique_ptr<cell>* __restrict__ cells)
{
	real_t position_difference[dims];

	real_t distance = velocity_helper<dims>::difference_and_distance(position + lhs * dims, position + rhs * dims,
																	 position_difference);

	const real_t adhesion_distance =
		relative_maximum_adhesion_distance[lhs] * radius[lhs] + relative_maximum_adhesion_distance[rhs] * radius[rhs];

	if (distance > adhesion_distance)
	{
		return;
	}

	cells[lhs]->neighbors().push_back(cells[rhs].get());
	cells[rhs]->neighbors().push_back(cells[lhs].get());

	// compute repulsion
	real_t repulsion;
	{
		const real_t repulsive_distance = radius[lhs] + radius[rhs];

		repulsion = 1 - distance / repulsive_distance;

		repulsion *= repulsion;

		repulsion *= sqrt(cell_cell_repulsion_strength[lhs] * cell_cell_repulsion_strength[rhs]);
	}

	simple_pressure[lhs] += repulsion * simple_pressure_coefficient;
	simple_pressure[rhs] += repulsion * simple_pressure_coefficient;

	// compute adhesion
	real_t adhesion;
	{
		adhesion = 1 - distance / adhesion_distance;

		adhesion *= adhesion;

		adhesion *= sqrt(cell_cell_adhesion_strength[lhs] * cell_cell_adhesion_strength[rhs]
						 * lhs_cell_adhesion_affinity * rhs_cell_adhesion_affinity);
	}

	real_t force = (repulsion - adhesion) / distance;

	velocity_helper<dims>::update_velocities(velocity + lhs * dims, velocity + rhs * dims, position_difference, force);
}

template <index_t dims>
void solve_internal(environment& e, index_t agents_count, index_t cell_def_count, real_t* __restrict__ velocity,
					real_t* __restrict__ simple_pressure, const real_t* __restrict__ position,
					const real_t* __restrict__ radius, const real_t* __restrict__ cell_cell_repulsion_strength,
					const real_t* __restrict__ cell_cell_adhesion_strength,
					const real_t* __restrict__ relative_maximum_adhesion_distance,
					const index_t* __restrict__ cell_definition_index,
					const real_t* __restrict__ cell_adhesion_affinities)
{
	auto cell_data = e.cells().cells().data();

	for (index_t i = 0; i < agents_count; i++)
	{
		common_solver::for_each_in_mech_neighborhood(
			e, common_solver::get_mesh_position(position + dims * i, e.mechanics_mesh), i, [=](index_t j) {
				const index_t i_cell_def_index = cell_definition_index[i];
				const index_t j_cell_def_index = cell_definition_index[j];

				solve_pair<dims>(i, j, velocity, simple_pressure, position, radius, cell_cell_repulsion_strength,
								 cell_cell_adhesion_strength, relative_maximum_adhesion_distance,
								 cell_adhesion_affinities[i * cell_def_count + j_cell_def_index],
								 cell_adhesion_affinities[j * cell_def_count + i_cell_def_index], cell_data);
			});
	}
}

void velocity_solver::solve(environment& e)
{
	auto& data = get_cell_data(e);

	clear_simple_pressure(data.simple_pressures.data(), data.agents_count);

	// clear neighbors
	for (index_t i = 0; i < data.agents_count; i++)
		e.cells().cells()[i]->neighbors().clear();

	if (e.m.mesh.dims == 1)
		solve_internal<1>(e, data.agents_count, e.cell_definitions_count, data.velocities.data(),
						  data.simple_pressures.data(), data.agent_data.positions.data(), data.geometries.radius.data(),
						  data.mechanics.cell_cell_repulsion_strength.data(),
						  data.mechanics.cell_cell_adhesion_strength.data(),
						  data.mechanics.relative_maximum_adhesion_distance.data(), data.cell_definition_indices.data(),
						  data.mechanics.cell_adhesion_affinities.data());
	else if (e.m.mesh.dims == 2)
		solve_internal<2>(e, data.agents_count, e.cell_definitions_count, data.velocities.data(),
						  data.simple_pressures.data(), data.agent_data.positions.data(), data.geometries.radius.data(),
						  data.mechanics.cell_cell_repulsion_strength.data(),
						  data.mechanics.cell_cell_adhesion_strength.data(),
						  data.mechanics.relative_maximum_adhesion_distance.data(), data.cell_definition_indices.data(),
						  data.mechanics.cell_adhesion_affinities.data());
	else if (e.m.mesh.dims == 3)
		solve_internal<3>(e, data.agents_count, e.cell_definitions_count, data.velocities.data(),
						  data.simple_pressures.data(), data.agent_data.positions.data(), data.geometries.radius.data(),
						  data.mechanics.cell_cell_repulsion_strength.data(),
						  data.mechanics.cell_cell_adhesion_strength.data(),
						  data.mechanics.relative_maximum_adhesion_distance.data(), data.cell_definition_indices.data(),
						  data.mechanics.cell_adhesion_affinities.data());
}
