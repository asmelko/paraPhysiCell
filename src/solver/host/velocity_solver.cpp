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
void solve_pair(index_t lhs, index_t rhs, real_t* __restrict__ velocity, real_t* __restrict__ simple_pressure,
				const real_t* __restrict__ position, const real_t* __restrict__ radius,
				const real_t* __restrict__ cell_cell_repulsion_strength,
				const real_t* __restrict__ cell_cell_adhesion_strength,
				const real_t* __restrict__ relative_maximum_adhesion_distance, real_t lhs_cell_adhesion_affinity,
				real_t rhs_cell_adhesion_affinity)
{
	real_t position_difference[dims];

	real_t distance = 0;

	for (index_t i = 0; i < dims; i++)
	{
		position_difference[i] = position[lhs * dims + i] - position[rhs * dims + i];
		distance += position_difference[i] * position_difference[i];
	}

	distance = std::sqrt(distance);

	// compute repulsion
	real_t repulsion;
	{
		const real_t repulsive_distance = radius[lhs] + radius[rhs];

		repulsion = 1 - distance / repulsive_distance;

		repulsion = repulsion > 0 ? repulsion : 0;

		repulsion *= repulsion;

		repulsion *= sqrt(cell_cell_repulsion_strength[lhs] * cell_cell_adhesion_strength[rhs]);
	}

	simple_pressure[lhs] += repulsion * simple_pressure_coefficient;
	simple_pressure[rhs] += repulsion * simple_pressure_coefficient;

	// compute adhesion
	real_t adhesion;
	{
		const real_t adhesion_distance = relative_maximum_adhesion_distance[lhs] * radius[lhs]
										 + relative_maximum_adhesion_distance[rhs] * radius[rhs];

		adhesion = 1 - distance / adhesion_distance;

		adhesion = adhesion > 0 ? adhesion : 0;

		adhesion *= adhesion;

		adhesion *= sqrt(cell_cell_adhesion_strength[lhs] * cell_cell_adhesion_strength[rhs]
						 * lhs_cell_adhesion_affinity * rhs_cell_adhesion_affinity);
	}

	real_t force = (repulsion - adhesion) / distance;

	for (index_t i = 0; i < dims; i++)
	{
		velocity[lhs * dims + i] += force * position_difference[i];
		velocity[rhs * dims + i] -= force * position_difference[i];
	}
}

template <index_t dims>
void solve_internal(index_t agents_count, index_t cell_def_count, real_t* __restrict__ velocity,
					real_t* __restrict__ simple_pressure, const real_t* __restrict__ position,
					const real_t* __restrict__ radius, const real_t* __restrict__ cell_cell_repulsion_strength,
					const real_t* __restrict__ cell_cell_adhesion_strength,
					const real_t* __restrict__ relative_maximum_adhesion_distance,
					const index_t* __restrict__ cell_definition_index,
					const real_t* __restrict__ cell_adhesion_affinities)
{
	for (index_t i = 0; i < agents_count; i++)
	{
		for (index_t j = i + 1; j < agents_count; j++)
		{
			if (cell_definition_index[i] == cell_definition_index[j])
			{
				const index_t i_cell_index = cell_definition_index[i];
				const index_t j_cell_index = cell_definition_index[j];

				solve_pair<dims>(i, j, velocity, simple_pressure, position, radius, cell_cell_repulsion_strength,
								 cell_cell_adhesion_strength, relative_maximum_adhesion_distance,
								 cell_adhesion_affinities[i * cell_def_count + j_cell_index],
								 cell_adhesion_affinities[j * cell_def_count + i_cell_index]);
			}
		}
	}
}

void velocity_solver::solve(cell_data& data)
{
	clear_simple_pressure(data.simple_pressure.data(), data.agents_count);

	if (data.m.mesh.dims == 1)
		solve_internal<1>(data.agents_count, data.mechanics.cell_definitions_count, data.velocities.data(),
						  data.simple_pressure.data(), data.agent_data.positions.data(), data.geometry.radius.data(),
						  data.mechanics.cell_cell_repulsion_strength.data(),
						  data.mechanics.cell_cell_adhesion_strength.data(),
						  data.mechanics.relative_maximum_adhesion_distance.data(), data.cell_definition_index.data(),
						  data.mechanics.cell_adhesion_affinities.data());
	else if (data.m.mesh.dims == 2)
		solve_internal<2>(data.agents_count, data.mechanics.cell_definitions_count, data.velocities.data(),
						  data.simple_pressure.data(), data.agent_data.positions.data(), data.geometry.radius.data(),
						  data.mechanics.cell_cell_repulsion_strength.data(),
						  data.mechanics.cell_cell_adhesion_strength.data(),
						  data.mechanics.relative_maximum_adhesion_distance.data(), data.cell_definition_index.data(),
						  data.mechanics.cell_adhesion_affinities.data());
	else if (data.m.mesh.dims == 3)
		solve_internal<3>(data.agents_count, data.mechanics.cell_definitions_count, data.velocities.data(),
						  data.simple_pressure.data(), data.agent_data.positions.data(), data.geometry.radius.data(),
						  data.mechanics.cell_cell_repulsion_strength.data(),
						  data.mechanics.cell_cell_adhesion_strength.data(),
						  data.mechanics.relative_maximum_adhesion_distance.data(), data.cell_definition_index.data(),
						  data.mechanics.cell_adhesion_affinities.data());
}
