#include "morse_position_solver.h"

#include <BioFVM/microenvironment.h>

#include "../../models/morse_position_model.h"
#include "../../random.h"
#include "solver_helper.h"

using namespace biofvm;
using namespace physicell;

template <index_t dims>
void solve_pair_intra(index_t lhs, index_t rhs, real_t* __restrict__ velocity, const real_t* __restrict__ position,
					  const real_t* __restrict__ scaling_factors, const real_t* __restrict__ equilibrium_distances,
					  const real_t* __restrict__ stiffnesses)
{
	real_t scaling_factor = scaling_factors[lhs];
	real_t equilibrium_distance = equilibrium_distances[lhs];
	real_t stiffness = stiffnesses[lhs];

	real_t potential_well_depth =
		(stiffness * equilibrium_distance * equilibrium_distance) / (8 * scaling_factor * scaling_factor);

	real_t position_difference[dims];

	const real_t distance = position_helper<dims>::difference_and_distance(position + lhs * dims, position + rhs * dims,
																		   position_difference);

	// if (distance > 2.5 * equilibrium_distance)
	// 	return;

	real_t exp_power = scaling_factor * (1 - (distance * distance) / (equilibrium_distance * equilibrium_distance));

	real_t force = (4 * scaling_factor * distance * potential_well_depth)
				   * (std::exp(exp_power) * std::exp(exp_power) - std::exp(exp_power))
				   / (equilibrium_distance * equilibrium_distance);

	position_helper<dims>::update_velocity(velocity + lhs * dims, position_difference, force / distance);

	// #pragma omp critical
	// 	{
	// 		std::cout << "intra Scaling factor: " << scaling_factor << " Equilibrium distance: " << equilibrium_distance
	// 				  << " Stiffness: " << stiffness << " Potential well depth: " << potential_well_depth
	// 				  << " Distance: " << distance << " Force: " << force << " Exp power:" << exp_power << std::endl;
	// 	}
}

template <index_t dims>
void solve_pair_inter(index_t lhs, index_t rhs, index_t cell_defs_count, real_t* __restrict__ velocity,
					  const real_t* __restrict__ position, const real_t* __restrict__ scaling_factors,
					  const real_t* __restrict__ equilibrium_distances, const real_t* __restrict__ stiffnesses,
					  const index_t* __restrict__ cell_definition_index)
{
	const auto lhs_index = cell_definition_index[lhs];
	const auto rhs_index = cell_definition_index[rhs];

	real_t scaling_factor = scaling_factors[cell_defs_count * lhs_index + rhs_index];
	real_t equilibrium_distance = equilibrium_distances[cell_defs_count * lhs_index + rhs_index];
	real_t stiffness = stiffnesses[cell_defs_count * lhs_index + rhs_index];

	real_t potential_well_depth =
		(stiffness * equilibrium_distance * equilibrium_distance) / (8 * scaling_factor * scaling_factor);

	real_t position_difference[dims];

	const real_t distance = position_helper<dims>::difference_and_distance(position + lhs * dims, position + rhs * dims,
																		   position_difference);

	real_t exp_power = scaling_factor * (1 - (distance * distance) / (equilibrium_distance * equilibrium_distance));

	real_t force = (4 * scaling_factor * distance * potential_well_depth)
				   * (std::exp(exp_power) * std::exp(exp_power) - std::exp(exp_power))
				   / (equilibrium_distance * equilibrium_distance);

	position_helper<dims>::update_velocity(velocity + lhs * dims, position_difference, force / distance);

	// #pragma omp critical
	// 	{
	// 		std::cout << "inter Scaling factor: " << scaling_factor << " Equilibrium distance: " << equilibrium_distance
	// 				  << " Stiffness: " << stiffness << " Potential well depth: " << potential_well_depth
	// 				  << " Distance: " << distance << " Force: " << force << " Exp power:" << exp_power << std::endl;
	// 	}
}

template <index_t dims>
void update_motility_single(index_t i, real_t time_step, real_t* __restrict__ velocity,
							const real_t* __restrict__ persistence_time, const std::uint8_t* __restrict__ is_motile,
							const real_t* __restrict__ migration_speed, const real_t* __restrict__ viscosity)
{
	if (is_motile[i] == 0)
		return;

	if (random::instance().uniform() < time_step / persistence_time[i])
	{
		real_t random_walk[dims];

		for (index_t d = 0; d < dims; d++)
			random_walk[d] = random::instance().normal(0, migration_speed[i] * viscosity[i]);

		position_helper<dims>::update_velocity(velocity + i * dims, random_walk, 1);
	}
}

template <index_t dims>
void update_motility_internal(index_t agents_count, real_t time_step, real_t* __restrict__ velocity,
							  const real_t* __restrict__ persistence_time, const std::uint8_t* __restrict__ is_motile,
							  const real_t* __restrict__ migration_speed, const real_t* __restrict__ viscosity)
{
#pragma omp for
	for (index_t i = 0; i < agents_count; i++)
	{
		update_motility_single<dims>(i, time_step, velocity, persistence_time, is_motile, migration_speed, viscosity);
	}
}

template <index_t dims>
void update_cell_forces_internal(cell_data& data, index_t cell_definitions_count, environment& e)
{
#pragma omp for
	for (index_t i = 0; i < data.agents_count; i++)
	{
		if (data.is_movable[i] == 0)
			continue;

		for (index_t j = 0; j < data.agents_count; j++)
		{
			if (i == j)
				continue;

			if (data.cell_residency[i] == data.cell_residency[j])
			{
				solve_pair_intra<dims>(i, j, data.velocities.data(), data.agent_data.positions.data(),
									   data.intra_scaling_factors.data(), data.intra_equilibrium_distances.data(),
									   data.intra_stiffnesses.data());
			}
			else
			{
				solve_pair_inter<dims>(i, j, cell_definitions_count, data.velocities.data(),
									   data.agent_data.positions.data(), e.inter_scaling_factors.data(),
									   e.inter_equilibrium_distances.data(), e.inter_stiffnesses.data(),
									   data.cell_definition_indices.data());
			}
		}
	}
}

void morse_position_solver::update_cell_forces(environment& e)
{
	auto& data = get_cell_data(e);

	update_cell_forces_internal<2>(data, e.cell_definitions_count, e);
}

void morse_position_solver::update_motility(environment& e)
{
	auto& data = get_cell_data(e);

	if (e.m.mesh.dims == 1)
		update_motility_internal<1>(data.agents_count, e.mechanics_time_step, data.velocities.data(),
									data.motilities.persistence_time.data(), data.motilities.is_motile.data(),
									data.motilities.migration_speed.data(), data.viscosities.data());
	else if (e.m.mesh.dims == 2)
		update_motility_internal<2>(data.agents_count, e.mechanics_time_step, data.velocities.data(),
									data.motilities.persistence_time.data(), data.motilities.is_motile.data(),
									data.motilities.migration_speed.data(), data.viscosities.data());
	else if (e.m.mesh.dims == 3)
		update_motility_internal<3>(data.agents_count, e.mechanics_time_step, data.velocities.data(),
									data.motilities.persistence_time.data(), data.motilities.is_motile.data(),
									data.motilities.migration_speed.data(), data.viscosities.data());
}

template <index_t dims>
void update_positions_internal(index_t agents_count, real_t time_step, real_t* __restrict__ position,
							   real_t* __restrict__ velocity, real_t* __restrict__ prev_velocity,
							   const real_t* __restrict__ viscosities, cell_data& data)
{
#pragma omp for
	for (index_t i = 0; i < agents_count; i++)
	{
		if (!data.is_movable[i])
			continue;

		point_t<real_t, 2> copy_vel;
		copy_vel[0] = velocity[i * dims + 0];
		copy_vel[1] = velocity[i * dims + 1];

		velocity[i * dims + 0] = 0;
		velocity[i * dims + 1] = 0;

		prev_velocity[i * dims + 0] = copy_vel[0];
		prev_velocity[i * dims + 1] = copy_vel[1];

		data.prev_velocities[i].push_back(copy_vel);

		if (data.prev_velocities[i].size() > 1)
			data.prev_velocities[i].erase(data.prev_velocities[i].begin());

		for (index_t d = 0; d < dims; d++)
		{
			if (data.prev_velocities[i].size() == 1)
			{
				position[i * dims + d] += (1 / viscosities[i]) * time_step * data.prev_velocities[i][0][d];
			}
			if (data.prev_velocities[i].size() == 2)
			{
				position[i * dims + d] +=
					(1 / viscosities[i]) * time_step
					* (-0.5 * data.prev_velocities[i][0][d] + 1.5 * data.prev_velocities[i][1][d]);
			}
			if (data.prev_velocities[i].size() == 3)
			{
				position[i * dims + d] +=
					(1 / viscosities[i]) * time_step
					* ((5. / 12.) * data.prev_velocities[i][0][d] + (-16. / 12) * data.prev_velocities[i][1][d]
					   + (23. / 12) * data.prev_velocities[i][2][d]);
			}
			if (data.prev_velocities[i].size() == 4)
			{
				position[i * dims + d] +=
					(1 / viscosities[i]) * time_step
					* ((-3. / 8.) * data.prev_velocities[i][0][d] + (37. / 24) * data.prev_velocities[i][1][d]
					   + (-59. / 24) * data.prev_velocities[i][2][d] + (55. / 24) * data.prev_velocities[i][3][d]);
			}
		}

		// #pragma omp critical
		// 		{
		// 			std::cout << "Velocity[0]: " << copy_vel[0] << " Velocity[1]: " << copy_vel[1]
		// 					  << " Viscosity: " << viscosities[i] << " position[0]: " << position[i * dims + 0]
		// 					  << " position[1]: " << position[i * dims + 1] << " velocity[0]/viscosity "
		// 					  << copy_vel[0] / (viscosities[i] / agents_count) << " velocity[1]/viscosity "
		// 					  << copy_vel[1] / (viscosities[i] / agents_count) << std::endl;
		// 		}
	}
}

void morse_position_solver::update_positions(environment& e)
{
	update_positions_internal<2>(get_cell_data(e).agents_count, e.mechanics_time_step,
								 get_cell_data(e).agent_data.positions.data(), get_cell_data(e).velocities.data(),
								 get_cell_data(e).previous_velocities.data(), get_cell_data(e).viscosities.data(),
								 get_cell_data(e));
}

void morse_position_model::update_cell_forces(environment& e) { morse_position_solver::update_cell_forces(e); }

void morse_position_model::update_motility(environment& e) { morse_position_solver::update_motility(e); }

void morse_position_model::update_positions(environment& e) { morse_position_solver::update_positions(e); }
