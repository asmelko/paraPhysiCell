#include "position_solver.h"

#include <cmath>
#include <cstdint>

#include <BioFVM/microenvironment.h>

#include "../../random.h"
#include "mesh.h"
#include "types.h"

using namespace biofvm;
using namespace physicell;

constexpr real_t simple_pressure_coefficient = 36.64504274775163; // 1 / (12 * (1 - sqrt(pi/(2*sqrt(3))))^2)
constexpr real_t zero_threshold = 1e-16;

void clear_simple_pressure(real_t* __restrict__ simple_pressure, index_t count)
{
	for (index_t i = 0; i < count; i++)
	{
		simple_pressure[i] = 0;
	}
}

constexpr void update_membrane_velocity(real_t position, real_t bounding_box, real_t sign, real_t radius,
										real_t repulsion_strength, real_t& velocity)
{
	const real_t dist = std::abs(bounding_box - position);

	real_t repulsion = 1 - dist / radius;
	repulsion = repulsion < 0 ? repulsion : 0;

	repulsion *= repulsion * repulsion_strength * sign;

	velocity += repulsion * dist;
}

template <index_t dims>
struct position_helper
{};

template <>
struct position_helper<1>
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

	static void random_walk(bool, real_t* __restrict__ walk)
	{
		real_t rand = random::instance().uniform();
		walk[0] = rand < 0.5 ? -1 : 1;
	}

	static constexpr void update_motility_vector(real_t* __restrict__ motility_vector, const real_t* __restrict__ walk,
												 const real_t* __restrict__ migration_bias_direction,
												 const real_t migration_bias)
	{
		motility_vector[0] = (1 - migration_bias) * walk[0] + migration_bias * migration_bias_direction[0];
	}

	static constexpr void normalize_and_scale(real_t* __restrict__ vector, real_t scale)
	{
		real_t length = std::abs(vector[0]);

		vector[0] = length > zero_threshold ? vector[0] * scale / length : 0;
	}

	static constexpr void update_membrane_velocities(real_t* __restrict__ velocity, const real_t* __restrict__ position,
													 const cartesian_mesh& mesh, const real_t radius,
													 const real_t repulsion_strength)
	{
		update_membrane_velocity(position[0], mesh.bounding_box_mins[0], 1, radius, repulsion_strength, velocity[0]);
		update_membrane_velocity(position[0], mesh.bounding_box_maxs[1], -1, radius, repulsion_strength, velocity[0]);
	}

	static constexpr void add(real_t* __restrict__ lhs, const real_t* __restrict__ rhs) { lhs[0] += rhs[0]; }
};

template <>
struct position_helper<2>
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

	static void random_walk(bool, real_t* __restrict__ walk)
	{
		real_t theta = random::instance().uniform(0, 2 * M_PI);
		walk[0] = std::cos(theta);
		walk[1] = std::sin(theta);
	}

	static constexpr void update_motility_vector(real_t* __restrict__ motility_vector, const real_t* __restrict__ walk,
												 const real_t* __restrict__ migration_bias_direction,
												 const real_t migration_bias)
	{
		motility_vector[0] = (1 - migration_bias) * walk[0] + migration_bias * migration_bias_direction[0];
		motility_vector[1] = (1 - migration_bias) * walk[1] + migration_bias * migration_bias_direction[1];
	}

	static constexpr void normalize_and_scale(real_t* __restrict__ vector, real_t scale)
	{
		real_t length = std::sqrt(vector[0] * vector[0] + vector[1] * vector[1]);

		vector[0] = length > zero_threshold ? vector[0] * scale / length : 0;
		vector[1] = length > zero_threshold ? vector[1] * scale / length : 0;
	}

	static constexpr void update_membrane_velocities(real_t* __restrict__ velocity, const real_t* __restrict__ position,
													 const cartesian_mesh& mesh, const real_t radius,
													 const real_t repulsion_strength)
	{
		update_membrane_velocity(position[0], mesh.bounding_box_mins[0], 1, radius, repulsion_strength, velocity[0]);
		update_membrane_velocity(position[0], mesh.bounding_box_maxs[0], -1, radius, repulsion_strength, velocity[0]);
		update_membrane_velocity(position[1], mesh.bounding_box_mins[1], 1, radius, repulsion_strength, velocity[1]);
		update_membrane_velocity(position[1], mesh.bounding_box_maxs[1], -1, radius, repulsion_strength, velocity[1]);
	}

	static constexpr void add(real_t* __restrict__ lhs, const real_t* __restrict__ rhs)
	{
		lhs[0] += rhs[0];
		lhs[1] += rhs[1];
	}
};

template <>
struct position_helper<3>
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

	static void random_walk(bool restrict_to_2d, real_t* __restrict__ walk)
	{
		if (restrict_to_2d)
		{
			position_helper<2>::random_walk(true, walk);
			walk[2] = 0;
		}
		else
		{
			const real_t theta = random::instance().uniform(0, 2 * M_PI);
			const real_t z = random::instance().uniform(-1, 1);
			const real_t r = std::sqrt(1 - z * z);

			walk[0] = std::cos(theta) * r;
			walk[1] = std::sin(theta) * r;
			walk[2] = z;
		}
	}

	static constexpr void update_motility_vector(real_t* __restrict__ motility_vector, const real_t* __restrict__ walk,
												 const real_t* __restrict__ migration_bias_direction,
												 const real_t migration_bias)
	{
		motility_vector[0] = (1 - migration_bias) * walk[0] + migration_bias * migration_bias_direction[0];
		motility_vector[1] = (1 - migration_bias) * walk[1] + migration_bias * migration_bias_direction[1];
		motility_vector[2] = (1 - migration_bias) * walk[2] + migration_bias * migration_bias_direction[2];
	}

	static constexpr void normalize_and_scale(real_t* __restrict__ vector, real_t scale)
	{
		real_t length = std::sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);

		vector[0] = length > zero_threshold ? vector[0] * scale / length : 0;
		vector[1] = length > zero_threshold ? vector[1] * scale / length : 0;
		vector[2] = length > zero_threshold ? vector[2] * scale / length : 0;
	}

	static constexpr void update_membrane_velocities(real_t* __restrict__ velocity, const real_t* __restrict__ position,
													 const cartesian_mesh& mesh, const real_t radius,
													 const real_t repulsion_strength)
	{
		update_membrane_velocity(position[0], mesh.bounding_box_mins[0], 1, radius, repulsion_strength, velocity[0]);
		update_membrane_velocity(position[0], mesh.bounding_box_maxs[0], -1, radius, repulsion_strength, velocity[0]);
		update_membrane_velocity(position[1], mesh.bounding_box_mins[1], 1, radius, repulsion_strength, velocity[1]);
		update_membrane_velocity(position[1], mesh.bounding_box_maxs[1], -1, radius, repulsion_strength, velocity[1]);
		update_membrane_velocity(position[2], mesh.bounding_box_mins[2], 1, radius, repulsion_strength, velocity[2]);
		update_membrane_velocity(position[2], mesh.bounding_box_maxs[2], -1, radius, repulsion_strength, velocity[2]);
	}

	static constexpr void add(real_t* __restrict__ lhs, const real_t* __restrict__ rhs)
	{
		lhs[0] += rhs[0];
		lhs[1] += rhs[1];
		lhs[2] += rhs[2];
	}
};

template <index_t dims>
void solve_pair(index_t lhs, index_t rhs, index_t cell_defs_count, real_t* __restrict__ velocity,
				real_t* __restrict__ simple_pressure, const real_t* __restrict__ position,
				const real_t* __restrict__ radius, const real_t* __restrict__ cell_cell_repulsion_strength,
				const real_t* __restrict__ cell_cell_adhesion_strength,
				const real_t* __restrict__ relative_maximum_adhesion_distance,
				const real_t* __restrict__ cell_adhesion_affinity, const index_t* __restrict__ cell_definition_index,
				std::unique_ptr<cell>* __restrict__ cells)
{
	real_t position_difference[dims];

	real_t distance = position_helper<dims>::difference_and_distance(position + lhs * dims, position + rhs * dims,
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

		const index_t lhs__cell_def_index = cell_definition_index[lhs];
		const index_t rhs__cell_def_index = cell_definition_index[rhs];

		adhesion *= sqrt(cell_cell_adhesion_strength[lhs] * cell_cell_adhesion_strength[rhs]
						 * cell_adhesion_affinity[lhs * cell_defs_count + rhs__cell_def_index]
						 * cell_adhesion_affinity[rhs * cell_defs_count + lhs__cell_def_index]);
	}

	real_t force = (repulsion - adhesion) / distance;

	position_helper<dims>::update_velocities(velocity + lhs * dims, velocity + rhs * dims, position_difference, force);
}

template <index_t dims>
void update_cell_forces_and_neighbors_single(
	environment& e, index_t i, index_t cell_def_count, real_t* __restrict__ velocity,
	real_t* __restrict__ simple_pressure, const real_t* __restrict__ position, const real_t* __restrict__ radius,
	const real_t* __restrict__ cell_cell_repulsion_strength, const real_t* __restrict__ cell_cell_adhesion_strength,
	const real_t* __restrict__ relative_maximum_adhesion_distance, const index_t* __restrict__ cell_definition_index,
	const real_t* __restrict__ cell_adhesion_affinities, std::unique_ptr<cell>* __restrict__ cell_data)
{
	common_solver::for_each_in_mech_neighborhood(
		e, common_solver::get_mesh_position(position + dims * i, e.mechanics_mesh), i, [=](index_t j) {
			solve_pair<dims>(i, j, cell_def_count, velocity, simple_pressure, position, radius,
							 cell_cell_repulsion_strength, cell_cell_adhesion_strength,
							 relative_maximum_adhesion_distance, cell_adhesion_affinities, cell_definition_index,
							 cell_data);
		});
}

template <index_t dims>
void update_motility_single(index_t i, real_t* __restrict__ motility_vector, real_t* __restrict__ velocity,
							const real_t* __restrict__ persistence_time, const real_t* __restrict__ migration_bias,
							const real_t* __restrict__ migration_bias_direction,
							const std::uint8_t* __restrict__ restrict_to_2d, const std::uint8_t* __restrict__ is_motile,
							const real_t* __restrict__ migration_speed, real_t time_step)
{
	if (is_motile[i] == 0)
		return;

	real_t rand = random::instance().uniform();

	if (rand >= time_step / persistence_time[i])
		return;

	real_t random_walk[dims];

	position_helper<dims>::random_walk(restrict_to_2d, random_walk);

	position_helper<dims>::update_motility_vector(motility_vector + i * dims, random_walk,
												  migration_bias_direction + i * dims, migration_bias[i]);

	position_helper<dims>::normalize_and_scale(motility_vector + i * dims, migration_speed[i]);

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
void update_cell_velocities_and_neighbors_internal(
	environment& e, index_t agents_count, index_t cell_def_count, real_t* __restrict__ velocity,
	real_t* __restrict__ simple_pressure, const real_t* __restrict__ position, const real_t* __restrict__ radius,
	const real_t* __restrict__ cell_cell_repulsion_strength, const real_t* __restrict__ cell_cell_adhesion_strength,
	const real_t* __restrict__ relative_maximum_adhesion_distance, const index_t* __restrict__ cell_definition_index,
	const real_t* __restrict__ cell_adhesion_affinities, const std::uint8_t* __restrict__ is_movable)
{
	auto cell_data = e.cells().cells().data();

	for (index_t i = 0; i < agents_count; i++)
	{
		if (is_movable[i] == 0)
			continue;

		update_cell_forces_and_neighbors_single<dims>(e, i, cell_def_count, velocity, simple_pressure, position, radius,
													  cell_cell_repulsion_strength, cell_cell_adhesion_strength,
													  relative_maximum_adhesion_distance, cell_definition_index,
													  cell_adhesion_affinities, cell_data);
	}
}

template <index_t dims>
void update_motility_internal(index_t agents_count, real_t* __restrict__ motility_vector, real_t* __restrict__ velocity,
							  const real_t* __restrict__ persistence_time, const real_t* __restrict__ migration_bias,
							  const real_t* __restrict__ migration_bias_direction,
							  const std::uint8_t* __restrict__ restrict_to_2d,
							  const std::uint8_t* __restrict__ is_motile, const real_t* __restrict__ migration_speed,
							  real_t time_step)
{
	for (index_t i = 0; i < agents_count; i++)
	{
		update_motility_single<dims>(i, motility_vector, velocity, persistence_time, migration_bias,
									 migration_bias_direction, restrict_to_2d, is_motile, migration_speed, time_step);
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
	for (index_t i = 0; i < agents_count; i++)
	{
		if (is_movable[i] == 0)
			continue;

		update_basement_membrane_interactions_single<dims>(i, velocity, position, radius, cell_BM_repulsion_strength,
														   mesh);
	}
}

void position_solver::update_cell_velocities_and_neighbors(environment& e)
{
	auto& data = get_cell_data(e);

	clear_simple_pressure(data.simple_pressures.data(), data.agents_count);

	// clear neighbors
	for (index_t i = 0; i < data.agents_count; i++)
		e.cells().cells()[i]->neighbors().clear();

	if (e.m.mesh.dims == 1)
		update_cell_velocities_and_neighbors_internal<1>(
			e, data.agents_count, e.cell_definitions_count, data.velocities.data(), data.simple_pressures.data(),
			data.agent_data.positions.data(), data.geometries.radius.data(),
			data.mechanics.cell_cell_repulsion_strength.data(), data.mechanics.cell_cell_adhesion_strength.data(),
			data.mechanics.relative_maximum_adhesion_distance.data(), data.cell_definition_indices.data(),
			data.mechanics.cell_adhesion_affinities.data(), data.is_movable.data());
	else if (e.m.mesh.dims == 2)
		update_cell_velocities_and_neighbors_internal<2>(
			e, data.agents_count, e.cell_definitions_count, data.velocities.data(), data.simple_pressures.data(),
			data.agent_data.positions.data(), data.geometries.radius.data(),
			data.mechanics.cell_cell_repulsion_strength.data(), data.mechanics.cell_cell_adhesion_strength.data(),
			data.mechanics.relative_maximum_adhesion_distance.data(), data.cell_definition_indices.data(),
			data.mechanics.cell_adhesion_affinities.data(), data.is_movable.data());
	else if (e.m.mesh.dims == 3)
		update_cell_velocities_and_neighbors_internal<3>(
			e, data.agents_count, e.cell_definitions_count, data.velocities.data(), data.simple_pressures.data(),
			data.agent_data.positions.data(), data.geometries.radius.data(),
			data.mechanics.cell_cell_repulsion_strength.data(), data.mechanics.cell_cell_adhesion_strength.data(),
			data.mechanics.relative_maximum_adhesion_distance.data(), data.cell_definition_indices.data(),
			data.mechanics.cell_adhesion_affinities.data(), data.is_movable.data());
}

void position_solver::update_motility(environment& e)
{
	auto& data = get_cell_data(e);

	if (e.m.mesh.dims == 1)
		update_motility_internal<1>(data.agents_count, data.motility.motility_vector.data(), data.velocities.data(),
									data.motility.persistence_time.data(), data.motility.migration_bias.data(),
									data.motility.migration_bias_direction.data(), data.motility.restrict_to_2d.data(),
									data.motility.is_motile.data(), data.motility.migration_speed.data(),
									e.mechanics_time_step);
	else if (e.m.mesh.dims == 2)
		update_motility_internal<2>(data.agents_count, data.motility.motility_vector.data(), data.velocities.data(),
									data.motility.persistence_time.data(), data.motility.migration_bias.data(),
									data.motility.migration_bias_direction.data(), data.motility.restrict_to_2d.data(),
									data.motility.is_motile.data(), data.motility.migration_speed.data(),
									e.mechanics_time_step);
	else if (e.m.mesh.dims == 3)
		update_motility_internal<3>(data.agents_count, data.motility.motility_vector.data(), data.velocities.data(),
									data.motility.persistence_time.data(), data.motility.migration_bias.data(),
									data.motility.migration_bias_direction.data(), data.motility.restrict_to_2d.data(),
									data.motility.is_motile.data(), data.motility.migration_speed.data(),
									e.mechanics_time_step);
}

void position_solver::update_basement_membrane_interactions(environment& e)
{
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
