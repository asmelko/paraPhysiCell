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
	real_t distance = std::abs(bounding_box - position);

	distance = std::max<real_t>(distance, 0.00001);

	real_t repulsion = 1 - distance / radius;
	repulsion = repulsion < 0 ? repulsion : 0;

	repulsion *= repulsion * repulsion_strength * sign;

	velocity += repulsion * distance;
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

	static constexpr void subtract(real_t* __restrict__ dst, const real_t* __restrict__ lhs,
								   const real_t* __restrict__ rhs)
	{
		dst[0] = lhs[0] - rhs[0];
	}

	static constexpr void update_position(real_t* __restrict__ position, const real_t* __restrict__ velocity,
										  const real_t factor, const real_t* __restrict__ previous_velocity,
										  const real_t previous_factor)
	{
		position[0] += velocity[0] * factor + previous_velocity[0] * previous_factor;
	}
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

	static constexpr void subtract(real_t* __restrict__ dst, const real_t* __restrict__ lhs,
								   const real_t* __restrict__ rhs)
	{
		dst[0] = lhs[0] - rhs[0];
		dst[1] = lhs[1] - rhs[1];
	}

	static constexpr void update_position(real_t* __restrict__ position, const real_t* __restrict__ velocity,
										  const real_t factor, const real_t* __restrict__ previous_velocity,
										  const real_t previous_factor)
	{
		position[0] += velocity[0] * factor + previous_velocity[0] * previous_factor;
		position[1] += velocity[1] * factor + previous_velocity[1] * previous_factor;
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

	static constexpr void subtract(real_t* __restrict__ dst, const real_t* __restrict__ lhs,
								   const real_t* __restrict__ rhs)
	{
		dst[0] = lhs[0] - rhs[0];
		dst[1] = lhs[1] - rhs[1];
		dst[2] = lhs[2] - rhs[2];
	}

	static constexpr void update_position(real_t* __restrict__ position, const real_t* __restrict__ velocity,
										  const real_t factor, const real_t* __restrict__ previous_velocity,
										  const real_t previous_factor)
	{
		position[0] += velocity[0] * factor + previous_velocity[0] * previous_factor;
		position[1] += velocity[1] * factor + previous_velocity[1] * previous_factor;
		position[2] += velocity[2] * factor + previous_velocity[2] * previous_factor;
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

	distance = std::max<real_t>(distance, 0.00001);

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

template <index_t dims>
void spring_contract_function(index_t agents_count, index_t cell_defs_count, real_t* __restrict__ velocity,
							  const index_t* __restrict__ cell_definition_index,
							  const real_t* __restrict__ attachment_elastic_constant,
							  const real_t* __restrict__ cell_adhesion_affinity, const real_t* __restrict__ position,
							  const std::uint8_t* __restrict__ is_movable,
							  const std::unique_ptr<cell>* __restrict__ cells)
{
	for (index_t this_cell_index = 0; this_cell_index < agents_count; this_cell_index++)
	{
		if (is_movable[this_cell_index] == 0)
			continue;

		for (index_t j = 0; j < (index_t)cells[this_cell_index]->spring_attached_cells().size(); j++)
		{
			const index_t other_cell_index = cells[this_cell_index]->spring_attached_cells()[j]->index();

			if (other_cell_index < this_cell_index)
				continue;

			const index_t this_cell_def_index = cell_definition_index[this_cell_index];
			const index_t other_cell_def_index = cell_definition_index[other_cell_index];

			const real_t adhesion =
				sqrt(attachment_elastic_constant[this_cell_index] * attachment_elastic_constant[other_cell_index]
					 * cell_adhesion_affinity[this_cell_index * cell_defs_count + other_cell_def_index]
					 * cell_adhesion_affinity[other_cell_index * cell_defs_count + this_cell_def_index]);

			real_t difference[dims];

			position_helper<dims>::subtract(difference, position + other_cell_index * dims,
											position + this_cell_index * dims);

			position_helper<dims>::update_velocities(velocity + this_cell_index * dims,
													 velocity + other_cell_index * dims, difference, adhesion);
		}
	}
}

void update_spring_attachments_internal(index_t agents_count, real_t time_step, index_t cell_defs_count,
										const real_t* __restrict__ detachment_rate,
										const real_t* __restrict__ attachment_rate,
										const real_t* __restrict__ cell_adhesion_affinities,
										const index_t* __restrict__ maximum_number_of_attachments,
										const index_t* __restrict__ cell_definition_index,
										const std::unique_ptr<cell>* __restrict__ cells)
{
	// detach cells from springs
	for (index_t i = 0; i < agents_count; i++)
	{
		auto this_cell = cells[i].get();

		for (index_t j = 0; j < (index_t)this_cell->spring_attached_cells().size(); j++)
		{
			if (random::instance().uniform() <= detachment_rate[i] * time_step)
			{
				auto other_cell = this_cell->spring_attached_cells()[j];

				this_cell->spring_attached_cells()[j] = this_cell->spring_attached_cells().back();
				this_cell->spring_attached_cells().pop_back();
				j--;

				auto it = std::find(other_cell->spring_attached_cells().begin(),
									other_cell->spring_attached_cells().end(), this_cell);

				*it = other_cell->spring_attached_cells().back();
				other_cell->spring_attached_cells().pop_back();
			}
		}
	}

	// attach cells to springs
	for (index_t this_cell_index = 0; this_cell_index < agents_count; this_cell_index++)
	{
		auto this_cell = cells[this_cell_index].get();

		for (std::size_t j = 0; j < this_cell->neighbors().size(); j++)
		{
			auto other_cell = this_cell->neighbors()[j];
			const auto other_cell_index = other_cell->index();

			const real_t affinity =
				cell_adhesion_affinities[this_cell_index * cell_defs_count + cell_definition_index[other_cell_index]];

			const real_t attachment_prob = attachment_rate[this_cell_index] * time_step * affinity;

			if (random::instance().uniform() <= attachment_prob)
			{
				if (std::find(this_cell->spring_attached_cells().begin(), this_cell->spring_attached_cells().end(),
							  other_cell)
					!= this_cell->spring_attached_cells().end())
					continue;

				if ((index_t)this_cell->spring_attached_cells().size()
					>= maximum_number_of_attachments[this_cell_index])
					break;

				if ((index_t)other_cell->spring_attached_cells().size()
					>= maximum_number_of_attachments[other_cell_index])
					continue;

				this_cell->spring_attached_cells().push_back(other_cell);
				other_cell->spring_attached_cells().push_back(this_cell);
			}
		}
	}
}

void position_solver::update_spring_attachments(environment& e)
{
	auto& data = get_cell_data(e);
	auto& cells = e.cells().cells();

	update_spring_attachments_internal(
		data.agents_count, e.mechanics_time_step, e.cell_definitions_count, data.mechanics.detachment_rate.data(),
		data.mechanics.attachment_rate.data(), data.mechanics.cell_adhesion_affinities.data(),
		data.mechanics.maximum_number_of_attachments.data(), data.cell_definition_indices.data(), cells.data());

	if (e.mechanics_mesh.dims == 1)
		spring_contract_function<1>(
			data.agents_count, e.cell_definitions_count, data.velocities.data(), data.cell_definition_indices.data(),
			data.mechanics.attachment_elastic_constant.data(), data.mechanics.cell_adhesion_affinities.data(),
			data.agent_data.positions.data(), data.is_movable.data(), cells.data());
	else if (e.mechanics_mesh.dims == 2)
		spring_contract_function<2>(
			data.agents_count, e.cell_definitions_count, data.velocities.data(), data.cell_definition_indices.data(),
			data.mechanics.attachment_elastic_constant.data(), data.mechanics.cell_adhesion_affinities.data(),
			data.agent_data.positions.data(), data.is_movable.data(), cells.data());
	else if (e.mechanics_mesh.dims == 3)
		spring_contract_function<3>(
			data.agents_count, e.cell_definitions_count, data.velocities.data(), data.cell_definition_indices.data(),
			data.mechanics.attachment_elastic_constant.data(), data.mechanics.cell_adhesion_affinities.data(),
			data.agent_data.positions.data(), data.is_movable.data(), cells.data());
}

template <index_t dims>
void update_positions_internal(index_t agents_count, index_t time_step, real_t* __restrict__ position,
							   real_t* __restrict__ velocity, real_t* __restrict__ previous_velocity,
							   const std::uint8_t* __restrict__ is_movable)
{
	for (index_t i = 0; i < agents_count; i++)
	{
		if (!is_movable[i])
			continue;

		position_helper<dims>::update_position(position + i * dims, velocity + i * dims, time_step * 1.5,
											   previous_velocity + i * dims, time_step * -0.5);

		for (index_t d = 0; d < dims; d++)
		{
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
