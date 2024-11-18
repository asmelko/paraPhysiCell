#include "position_solver.h"

#include <BioFVM/microenvironment.h>

#include "../../random.h"
#include "solver_helper.h"

using namespace biofvm;
using namespace physicell;

constexpr real_t simple_pressure_coefficient = 36.64504274775163; // 1 / (12 * (1 - sqrt(pi/(2*sqrt(3))))^2)

void clear_simple_pressure(real_t* __restrict__ simple_pressure, index_t count)
{
#pragma omp for
	for (index_t i = 0; i < count; i++)
	{
		simple_pressure[i] = 0;
	}
}

template <index_t dims>
void solve_pair(index_t lhs, index_t rhs, index_t cell_defs_count, real_t* __restrict__ velocity,
				real_t* __restrict__ simple_pressure, const real_t* __restrict__ position,
				const real_t* __restrict__ radius, const real_t* __restrict__ cell_cell_repulsion_strength,
				const real_t* __restrict__ cell_cell_adhesion_strength,
				const real_t* __restrict__ relative_maximum_adhesion_distance,
				const real_t* __restrict__ cell_adhesion_affinity, const index_t* __restrict__ cell_definition_index)
{
	real_t position_difference[dims];

	const real_t distance = std::max<real_t>(position_helper<dims>::difference_and_distance(
												 position + lhs * dims, position + rhs * dims, position_difference),
											 0.00001);

	// compute repulsion
	real_t repulsion;
	{
		const real_t repulsive_distance = radius[lhs] + radius[rhs];

		repulsion = 1 - distance / repulsive_distance;

		repulsion = repulsion < 0 ? 0 : repulsion;

		repulsion *= repulsion;

		// update simple pressure
		simple_pressure[lhs] += repulsion * simple_pressure_coefficient;
		simple_pressure[rhs] += repulsion * simple_pressure_coefficient;

		repulsion *= std::sqrt(cell_cell_repulsion_strength[lhs] * cell_cell_repulsion_strength[rhs]);
	}

	// compute adhesion
	real_t adhesion;
	{
		const real_t adhesion_distance = relative_maximum_adhesion_distance[lhs] * radius[lhs]
										 + relative_maximum_adhesion_distance[rhs] * radius[rhs];

		adhesion = 1 - distance / adhesion_distance;

		adhesion *= adhesion;

		const index_t lhs_cell_def_index = cell_definition_index[lhs];
		const index_t rhs_cell_def_index = cell_definition_index[rhs];

		adhesion *= std::sqrt(cell_cell_adhesion_strength[lhs] * cell_cell_adhesion_strength[rhs]
							  * cell_adhesion_affinity[lhs * cell_defs_count + rhs_cell_def_index]
							  * cell_adhesion_affinity[rhs * cell_defs_count + lhs_cell_def_index]);
	}

	real_t force = (repulsion - adhesion) / distance;

	position_helper<dims>::update_velocity(velocity + lhs * dims, position_difference, force);
}

template <index_t dims>
void update_cell_neighbors_single(environment& e, index_t i, const real_t* __restrict__ position,
								  const real_t* __restrict__ radius,
								  const real_t* __restrict__ relative_maximum_adhesion_distance,
								  std::vector<index_t>* __restrict__ neighbors)
{
	common_solver::for_each_in_mech_neighborhood(
		e, common_solver::get_mesh_position(position + dims * i, e.mechanics_mesh), i, [=](index_t j) {
			const real_t adhesion_distance =
				relative_maximum_adhesion_distance[i] * radius[i] + relative_maximum_adhesion_distance[j] * radius[j];

			const real_t distance = position_helper<dims>::distance(position + i * dims, position + j * dims);

			if (distance <= adhesion_distance)
			{
				neighbors[i].push_back(j);
			}
		});
}

template <index_t dims>
void update_cell_forces_single(
	index_t i, index_t cell_def_count, real_t* __restrict__ velocity, real_t* __restrict__ simple_pressure,
	const real_t* __restrict__ position, const real_t* __restrict__ radius,
	const real_t* __restrict__ cell_cell_repulsion_strength, const real_t* __restrict__ cell_cell_adhesion_strength,
	const real_t* __restrict__ relative_maximum_adhesion_distance, const index_t* __restrict__ cell_definition_index,
	const real_t* __restrict__ cell_adhesion_affinities, std::vector<index_t>* __restrict__ neighbors)
{
	for (const index_t j : neighbors[i])
	{
		solve_pair<dims>(i, j, cell_def_count, velocity, simple_pressure, position, radius,
						 cell_cell_repulsion_strength, cell_cell_adhesion_strength, relative_maximum_adhesion_distance,
						 cell_adhesion_affinities, cell_definition_index);
	}
}

template <index_t dims>
void update_motility_single(index_t i, real_t time_step, real_t* __restrict__ motility_vector,
							real_t* __restrict__ velocity, const real_t* __restrict__ persistence_time,
							const real_t* __restrict__ migration_bias, real_t* __restrict__ migration_bias_direction,
							const std::uint8_t* __restrict__ restrict_to_2d, const std::uint8_t* __restrict__ is_motile,
							const real_t* __restrict__ migration_speed,
							const motility_data::direction_update_func* __restrict__ update_migration_bias_direction_f,
							cell_container& cells)
{
	if (is_motile[i] == 0)
		return;

	if (random::instance().uniform() < time_step / persistence_time[i])
	{
		real_t random_walk[dims];

		position_helper<dims>::random_walk(restrict_to_2d, random_walk);

		if (update_migration_bias_direction_f[i] != nullptr)
		{
			update_migration_bias_direction_f[i](*cells.agents()[i]);
		}

		position_helper<dims>::update_motility_vector(motility_vector + i * dims, random_walk,
													  migration_bias_direction + i * dims, migration_bias[i]);

		position_helper<dims>::normalize_and_scale(motility_vector + i * dims, migration_speed[i]);
	}

	position_helper<dims>::add(velocity + i * dims, motility_vector + i * dims);
}

template <index_t dims>
void update_basement_membrane_interactions_single(index_t i, real_t* __restrict__ velocity,
												  const real_t* __restrict__ position,
												  const real_t* __restrict__ radius,
												  const real_t* __restrict__ cell_BM_repulsion_strength,
												  const cartesian_mesh& mesh)
{
	position_helper<dims>::update_membrane_velocities(velocity + i * dims, position + i * dims, mesh, radius[i],
													  cell_BM_repulsion_strength[i]);
}

template <index_t dims>
void update_cell_forces_internal(
	index_t agents_count, index_t cell_def_count, real_t* __restrict__ velocity, real_t* __restrict__ simple_pressure,
	const real_t* __restrict__ position, const real_t* __restrict__ radius,
	const real_t* __restrict__ cell_cell_repulsion_strength, const real_t* __restrict__ cell_cell_adhesion_strength,
	const real_t* __restrict__ relative_maximum_adhesion_distance, const index_t* __restrict__ cell_definition_index,
	const real_t* __restrict__ cell_adhesion_affinities, const std::uint8_t* __restrict__ is_movable,
	std::vector<index_t>* __restrict__ neighbors)
{
#pragma omp for
	for (index_t i = 0; i < agents_count; i++)
	{
		if (is_movable[i] == 0)
			continue;

		update_cell_forces_single<dims>(i, cell_def_count, velocity, simple_pressure, position, radius,
										cell_cell_repulsion_strength, cell_cell_adhesion_strength,
										relative_maximum_adhesion_distance, cell_definition_index,
										cell_adhesion_affinities, neighbors);
	}
}

template <index_t dims>
void update_cell_neighbors_internal(environment& e, index_t agents_count, const real_t* __restrict__ position,
									const real_t* __restrict__ radius,
									const real_t* __restrict__ relative_maximum_adhesion_distance,
									const std::uint8_t* __restrict__ is_movable,
									std::vector<index_t>* __restrict__ neighbors)
{
#pragma omp for
	for (index_t i = 0; i < agents_count; i++)
	{
		if (is_movable[i] == 0)
			continue;

		update_cell_neighbors_single<dims>(e, i, position, radius, relative_maximum_adhesion_distance, neighbors);
	}
}

template <index_t dims>
void update_motility_internal(
	index_t agents_count, real_t time_step, real_t* __restrict__ motility_vector, real_t* __restrict__ velocity,
	const real_t* __restrict__ persistence_time, const real_t* __restrict__ migration_bias,
	real_t* __restrict__ migration_bias_direction, const std::uint8_t* __restrict__ restrict_to_2d,
	const std::uint8_t* __restrict__ is_motile, const real_t* __restrict__ migration_speed,
	const motility_data::direction_update_func* __restrict__ update_migration_bias_direction_f, cell_container& cells)
{
#pragma omp for
	for (index_t i = 0; i < agents_count; i++)
	{
		update_motility_single<dims>(i, time_step, motility_vector, velocity, persistence_time, migration_bias,
									 migration_bias_direction, restrict_to_2d, is_motile, migration_speed,
									 update_migration_bias_direction_f, cells);
	}
}

template <index_t dims>
void update_basement_membrane_interactions_internal(index_t agents_count, real_t* __restrict__ velocity,
													const real_t* __restrict__ position,
													const real_t* __restrict__ radius,
													const real_t* __restrict__ cell_BM_repulsion_strength,
													const std::uint8_t* __restrict__ is_movable,
													const cartesian_mesh& mesh)
{
#pragma omp for
	for (index_t i = 0; i < agents_count; i++)
	{
		if (is_movable[i] == 0)
			continue;

		update_basement_membrane_interactions_single<dims>(i, velocity, position, radius, cell_BM_repulsion_strength,
														   mesh);
	}
}

void position_solver::update_cell_forces(environment& e)
{
	auto& data = get_cell_data(e);

	clear_simple_pressure(data.states.simple_pressure.data(), data.agents_count);

	if (e.m.mesh.dims == 1)
		update_cell_forces_internal<1>(
			data.agents_count, e.cell_definitions_count, data.velocities.data(), data.states.simple_pressure.data(),
			data.agent_data.positions.data(), data.geometries.radius.data(),
			data.mechanics.cell_cell_repulsion_strength.data(), data.mechanics.cell_cell_adhesion_strength.data(),
			data.mechanics.relative_maximum_adhesion_distance.data(), data.cell_definition_indices.data(),
			data.mechanics.cell_adhesion_affinities.data(), data.is_movable.data(), data.states.neighbors.data());
	else if (e.m.mesh.dims == 2)
		update_cell_forces_internal<2>(
			data.agents_count, e.cell_definitions_count, data.velocities.data(), data.states.simple_pressure.data(),
			data.agent_data.positions.data(), data.geometries.radius.data(),
			data.mechanics.cell_cell_repulsion_strength.data(), data.mechanics.cell_cell_adhesion_strength.data(),
			data.mechanics.relative_maximum_adhesion_distance.data(), data.cell_definition_indices.data(),
			data.mechanics.cell_adhesion_affinities.data(), data.is_movable.data(), data.states.neighbors.data());
	else if (e.m.mesh.dims == 3)
		update_cell_forces_internal<3>(
			data.agents_count, e.cell_definitions_count, data.velocities.data(), data.states.simple_pressure.data(),
			data.agent_data.positions.data(), data.geometries.radius.data(),
			data.mechanics.cell_cell_repulsion_strength.data(), data.mechanics.cell_cell_adhesion_strength.data(),
			data.mechanics.relative_maximum_adhesion_distance.data(), data.cell_definition_indices.data(),
			data.mechanics.cell_adhesion_affinities.data(), data.is_movable.data(), data.states.neighbors.data());
}

void position_solver::update_cell_neighbors(environment& e)
{
	auto& data = get_cell_data(e);

	// clear neighbors
#pragma omp for
	for (index_t i = 0; i < data.agents_count; i++)
		data.states.neighbors[i].clear();

	if (e.m.mesh.dims == 1)
		update_cell_neighbors_internal<1>(e, data.agents_count, data.agent_data.positions.data(),
										  data.geometries.radius.data(),
										  data.mechanics.relative_maximum_adhesion_distance.data(),
										  data.is_movable.data(), data.states.neighbors.data());
	else if (e.m.mesh.dims == 2)
		update_cell_neighbors_internal<2>(e, data.agents_count, data.agent_data.positions.data(),
										  data.geometries.radius.data(),
										  data.mechanics.relative_maximum_adhesion_distance.data(),
										  data.is_movable.data(), data.states.neighbors.data());
	else if (e.m.mesh.dims == 3)
		update_cell_neighbors_internal<3>(e, data.agents_count, data.agent_data.positions.data(),
										  data.geometries.radius.data(),
										  data.mechanics.relative_maximum_adhesion_distance.data(),
										  data.is_movable.data(), data.states.neighbors.data());
}

void position_solver::update_motility(environment& e)
{
	auto& data = get_cell_data(e);
	auto& cells = e.get_container();

	if (e.m.mesh.dims == 1)
		update_motility_internal<1>(
			data.agents_count, e.mechanics_time_step, data.motilities.motility_vector.data(), data.velocities.data(),
			data.motilities.persistence_time.data(), data.motilities.migration_bias.data(),
			data.motilities.migration_bias_direction.data(), data.motilities.restrict_to_2d.data(),
			data.motilities.is_motile.data(), data.motilities.migration_speed.data(),
			data.motilities.update_migration_bias_direction.data(), cells);
	else if (e.m.mesh.dims == 2)
		update_motility_internal<2>(
			data.agents_count, e.mechanics_time_step, data.motilities.motility_vector.data(), data.velocities.data(),
			data.motilities.persistence_time.data(), data.motilities.migration_bias.data(),
			data.motilities.migration_bias_direction.data(), data.motilities.restrict_to_2d.data(),
			data.motilities.is_motile.data(), data.motilities.migration_speed.data(),
			data.motilities.update_migration_bias_direction.data(), cells);
	else if (e.m.mesh.dims == 3)
		update_motility_internal<3>(
			data.agents_count, e.mechanics_time_step, data.motilities.motility_vector.data(), data.velocities.data(),
			data.motilities.persistence_time.data(), data.motilities.migration_bias.data(),
			data.motilities.migration_bias_direction.data(), data.motilities.restrict_to_2d.data(),
			data.motilities.is_motile.data(), data.motilities.migration_speed.data(),
			data.motilities.update_migration_bias_direction.data(), cells);
}

void position_solver::update_basement_membrane_interactions(environment& e)
{
	if (!e.virtual_wall_at_domain_edges)
		return;

	auto& data = get_cell_data(e);

	if (e.m.mesh.dims == 1)
		update_basement_membrane_interactions_internal<1>(
			data.agents_count, data.velocities.data(), data.agent_data.positions.data(), data.geometries.radius.data(),
			data.mechanics.cell_BM_repulsion_strength.data(), data.is_movable.data(), e.m.mesh);
	else if (e.m.mesh.dims == 2)
		update_basement_membrane_interactions_internal<2>(
			data.agents_count, data.velocities.data(), data.agent_data.positions.data(), data.geometries.radius.data(),
			data.mechanics.cell_BM_repulsion_strength.data(), data.is_movable.data(), e.m.mesh);
	else if (e.m.mesh.dims == 3)
		update_basement_membrane_interactions_internal<3>(
			data.agents_count, data.velocities.data(), data.agent_data.positions.data(), data.geometries.radius.data(),
			data.mechanics.cell_BM_repulsion_strength.data(), data.is_movable.data(), e.m.mesh);
}

template <index_t dims>
void spring_contract_function(index_t agents_count, index_t cell_defs_count, real_t* __restrict__ velocity,
							  const index_t* __restrict__ cell_definition_index,
							  const real_t* __restrict__ attachment_elastic_constant,
							  const real_t* __restrict__ cell_adhesion_affinity, const real_t* __restrict__ position,
							  const std::uint8_t* __restrict__ is_movable, std::vector<index_t>* __restrict__ springs)
{
#pragma omp for
	for (index_t this_cell_index = 0; this_cell_index < agents_count; this_cell_index++)
	{
		if (is_movable[this_cell_index] == 0)
			continue;

		for (std::size_t j = 0; j < springs[this_cell_index].size(); j++)
		{
			const index_t other_cell_index = springs[this_cell_index][j];

			const index_t this_cell_def_index = cell_definition_index[this_cell_index];
			const index_t other_cell_def_index = cell_definition_index[other_cell_index];

			const real_t adhesion =
				sqrt(attachment_elastic_constant[this_cell_index] * attachment_elastic_constant[other_cell_index]
					 * cell_adhesion_affinity[this_cell_index * cell_defs_count + other_cell_def_index]
					 * cell_adhesion_affinity[other_cell_index * cell_defs_count + this_cell_def_index]);

			real_t difference[dims];

			position_helper<dims>::subtract(difference, position + other_cell_index * dims,
											position + this_cell_index * dims);

			position_helper<dims>::update_velocity(velocity + this_cell_index * dims, difference, adhesion);
		}
	}
}

constexpr index_t erased_spring = -1;

void update_spring_attachments_internal(
	index_t agents_count, real_t time_step, index_t cell_defs_count, const real_t* __restrict__ detachment_rate,
	const real_t* __restrict__ attachment_rate, const real_t* __restrict__ cell_adhesion_affinities,
	const index_t* __restrict__ maximum_number_of_attachments, const index_t* __restrict__ cell_definition_index,
	const std::vector<index_t>* __restrict__ neighbors, std::vector<index_t>* __restrict__ springs)
{
// mark springs for detachment
#pragma omp for
	for (index_t this_cell_index = 0; this_cell_index < agents_count; this_cell_index++)
	{
		for (index_t j = 0; j < (index_t)springs[this_cell_index].size(); j++)
		{
			if (random::instance().uniform() <= detachment_rate[this_cell_index] * time_step)
			{
#pragma omp critical
				{
					const index_t other_cell_index = springs[this_cell_index][j];

					if (other_cell_index != erased_spring)
					{
						springs[this_cell_index][j] = erased_spring;

						*std::find(springs[other_cell_index].begin(), springs[other_cell_index].end(),
								   this_cell_index) = erased_spring;
					}
				}
			}
		}
	}

// remove marked springs
#pragma omp for
	for (index_t this_cell_index = 0; this_cell_index < agents_count; this_cell_index++)
	{
		auto it = std::remove(springs[this_cell_index].begin(), springs[this_cell_index].end(), erased_spring);

		springs[this_cell_index].erase(it, springs[this_cell_index].end());
	}

// attach cells to springs
#pragma omp for
	for (index_t this_cell_index = 0; this_cell_index < agents_count; this_cell_index++)
	{
		for (std::size_t j = 0; j < neighbors[this_cell_index].size(); j++)
		{
			const index_t other_cell_index = neighbors[this_cell_index][j];

			if (other_cell_index < this_cell_index)
				continue;

			const real_t affinity_l =
				cell_adhesion_affinities[this_cell_index * cell_defs_count + cell_definition_index[other_cell_index]];

			const real_t attachment_prob_l = attachment_rate[this_cell_index] * time_step * affinity_l;

			const real_t affinity_r =
				cell_adhesion_affinities[other_cell_index * cell_defs_count + cell_definition_index[this_cell_index]];

			const real_t attachment_prob_r = attachment_rate[other_cell_index] * time_step * affinity_r;

			if (random::instance().uniform() <= attachment_prob_l || random::instance().uniform() <= attachment_prob_r)
			{
#pragma omp critical
				{
					if ((index_t)springs[this_cell_index].size() < maximum_number_of_attachments[this_cell_index]
						&& (index_t)springs[other_cell_index].size() < maximum_number_of_attachments[other_cell_index])
					{
						springs[this_cell_index].push_back(other_cell_index);
						springs[other_cell_index].push_back(this_cell_index);
					}
				}
			}
		}
	}
}

void position_solver::update_spring_attachments(environment& e)
{
	if (!e.automated_spring_adhesion)
		return;

	auto& data = get_cell_data(e);

	update_spring_attachments_internal(
		data.agents_count, e.mechanics_time_step, e.cell_definitions_count, data.mechanics.detachment_rate.data(),
		data.mechanics.attachment_rate.data(), data.mechanics.cell_adhesion_affinities.data(),
		data.mechanics.maximum_number_of_attachments.data(), data.cell_definition_indices.data(),
		data.states.neighbors.data(), data.states.springs.data());

	if (e.mechanics_mesh.dims == 1)
		spring_contract_function<1>(
			data.agents_count, e.cell_definitions_count, data.velocities.data(), data.cell_definition_indices.data(),
			data.mechanics.attachment_elastic_constant.data(), data.mechanics.cell_adhesion_affinities.data(),
			data.agent_data.positions.data(), data.is_movable.data(), data.states.springs.data());
	else if (e.mechanics_mesh.dims == 2)
		spring_contract_function<2>(
			data.agents_count, e.cell_definitions_count, data.velocities.data(), data.cell_definition_indices.data(),
			data.mechanics.attachment_elastic_constant.data(), data.mechanics.cell_adhesion_affinities.data(),
			data.agent_data.positions.data(), data.is_movable.data(), data.states.springs.data());
	else if (e.mechanics_mesh.dims == 3)
		spring_contract_function<3>(
			data.agents_count, e.cell_definitions_count, data.velocities.data(), data.cell_definition_indices.data(),
			data.mechanics.attachment_elastic_constant.data(), data.mechanics.cell_adhesion_affinities.data(),
			data.agent_data.positions.data(), data.is_movable.data(), data.states.springs.data());
}

template <index_t dims>
void update_positions_internal(index_t agents_count, real_t time_step, real_t* __restrict__ position,
							   real_t* __restrict__ velocity, real_t* __restrict__ previous_velocity,
							   const std::uint8_t* __restrict__ is_movable)
{
#pragma omp for
	for (index_t i = 0; i < agents_count; i++)
	{
		if (!is_movable[i])
			continue;

		const real_t factor = time_step * 1.5;
		const real_t previous_factor = time_step * -0.5;

		for (index_t d = 0; d < dims; d++)
		{
			position[i * dims + d] +=
				velocity[i * dims + d] * factor + previous_velocity[i * dims + d] * previous_factor;

			previous_velocity[i * dims + d] = velocity[i * dims + d];
			velocity[i * dims + d] = 0;
		}
	}
}

void position_solver::update_positions(environment& e)
{
	if (e.mechanics_mesh.dims == 1)
		update_positions_internal<1>(get_cell_data(e).agents_count, e.mechanics_time_step,
									 get_cell_data(e).agent_data.positions.data(), get_cell_data(e).velocities.data(),
									 get_cell_data(e).previous_velocities.data(), get_cell_data(e).is_movable.data());
	else if (e.mechanics_mesh.dims == 2)
		update_positions_internal<2>(get_cell_data(e).agents_count, e.mechanics_time_step,
									 get_cell_data(e).agent_data.positions.data(), get_cell_data(e).velocities.data(),
									 get_cell_data(e).previous_velocities.data(), get_cell_data(e).is_movable.data());
	else if (e.mechanics_mesh.dims == 3)
		update_positions_internal<3>(get_cell_data(e).agents_count, e.mechanics_time_step,
									 get_cell_data(e).agent_data.positions.data(), get_cell_data(e).velocities.data(),
									 get_cell_data(e).previous_velocities.data(), get_cell_data(e).is_movable.data());
}
