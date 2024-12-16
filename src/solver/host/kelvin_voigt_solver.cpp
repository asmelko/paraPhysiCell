#include "kelvin_voigt_solver.h"

#include <BioFVM/microenvironment.h>

#include "../../models/kelvin_voigt_model.h"
#include "../../random.h"
#include "solver_helper.h"

using namespace biofvm;
using namespace physicell;

template <index_t dims>
void solve_pair_intra(index_t lhs, index_t rhs, real_t rest_length, real_t time_step, real_t* __restrict__ velocity,
					  const real_t* __restrict__ position, const real_t* __restrict__ spring_constants,
					  const real_t* __restrict__ dissipation_rates, const real_t* __restrict__ prev_velocity)
{
	real_t spring_constant = spring_constants[lhs];

	real_t position_difference[dims];

	const real_t distance = position_helper<dims>::difference_and_distance(position + rhs * dims, position + lhs * dims,
																		   position_difference);

	for (index_t i = 0; i < dims; i++)
		position_difference[i] = position_difference[i] / distance;

	real_t force = spring_constant * (distance - rest_length);

	real_t velocity_difference[dims];

	position_helper<dims>::subtract(velocity_difference, prev_velocity + rhs * dims, prev_velocity + lhs * dims);

	real_t damp_force =
		dissipation_rates[lhs] * time_step
		* (position_difference[0] * velocity_difference[0] + position_difference[1] * velocity_difference[1]);

	position_helper<dims>::update_velocity(velocity + lhs * dims, position_difference, force + damp_force);

	// #pragma omp critical
	// 	{
	// 		std::cout << "intra spring constant: " << spring_constant << " Equilibrium distance: " << rest_length
	// 				  << " Dissipation rate: " << dissipation_rates[lhs] << " Distance: " << distance << " Force: " <<
	// force
	// 				  << " Damp force: " << damp_force << " lhs prev_velocity: [" << prev_velocity[lhs * dims + 0] << ",
	// "
	// 				  << prev_velocity[lhs * dims + 1] << "] rhs prev_velocity: [" << prev_velocity[rhs * dims + 0] <<
	// ", "
	// 				  << prev_velocity[rhs * dims + 1] << "]" << std::endl;
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
void update_cell_forces_internal(cell_data& data, index_t cell_definitions_count, environment& e)
{
#pragma omp for
	for (index_t i = 0; i < data.agents_count; i++)
	{
		if (data.is_movable[i] == 0)
			continue;

		for (index_t nei_i = 0; nei_i < (index_t)data.states.neighbors[i].size(); nei_i++)
		{
			solve_pair_intra<dims>(i, data.states.neighbors[i][nei_i], data.rest_lengths[i][nei_i],
								   e.mechanics_time_step, data.velocities.data(), data.agent_data.positions.data(),
								   data.spring_constants.data(), data.dissipation_rates.data(),
								   data.previous_velocities.data());
		}

		// experimental area preservation
		// for (index_t nei_i = 0; nei_i < (index_t)data.states.neighbors[i].size(); nei_i++)
		// {
		// 	for (index_t nei_j = nei_i + 1; nei_j < (index_t)data.states.neighbors[i].size(); nei_j++)
		// 	{
		// 		auto nei_i_idx = data.states.neighbors[i][nei_i];
		// 		auto nei_j_idx = data.states.neighbors[i][nei_j];

		// 		const real_t* pos = data.agent_data.positions.data() + i * dims;
		// 		const real_t* pos_nei_i = data.agent_data.positions.data() + nei_i_idx * dims;
		// 		const real_t* pos_nei_j = data.agent_data.positions.data() + nei_j_idx * dims;

		// 		if (std::find(data.states.neighbors[nei_i_idx].begin(), data.states.neighbors[nei_i_idx].end(),
		// 					  nei_j_idx)
		// 			== data.states.neighbors[nei_i_idx].end())
		// 			continue;

		// 		const real_t area = pos[0] * (pos_nei_i[1] - pos_nei_j[1]) + pos_nei_i[0] * (pos_nei_j[1] - pos[1])
		// 							+ pos_nei_j[0] * (pos[1] - pos_nei_i[1]) + 0.00000001;

		// 		// #pragma omp critical
		// 		// 				{
		// 		// 					std::cout << "Area of triangle " << i << " " << nei_i_idx << " " << nei_j_idx << "
		// 		// is " << area
		// 		// 							  << std::endl;
		// 		// 				}

		// 		data.velocities[i * dims + 0] += 0.00001 * (pos_nei_i[1] - pos_nei_j[1]) * std::pow(area, 1);
		// 		data.velocities[i * dims + 1] += 0.00001 * (pos_nei_j[0] - pos_nei_i[0]) * std::pow(area, 1);
		// 	}
		// }

		for (index_t j = 0; j < data.agents_count; j++)
		{
			if (i == j)
				continue;

			if (data.cell_residency[i] != data.cell_residency[j])
			{
				solve_pair_inter<dims>(i, j, cell_definitions_count, data.velocities.data(),
									   data.agent_data.positions.data(), e.inter_scaling_factors.data(),
									   e.inter_equilibrium_distances.data(), e.inter_stiffnesses.data(),
									   data.cell_definition_indices.data());
			}
		}
	}
}

bool cells_updated = false;

template <index_t dims>
void update_cell_neighbors_internal(cell_data& data)
{
	if (cells_updated)
	{
		return;
	}

	cells_updated = true;

	for (index_t i = 0; i < data.agents_count; i++)
	{
		for (index_t j = 0; j < data.agents_count; j++)
		{
			if (i == j || data.cell_residency[i] != data.cell_residency[j])
				continue;

			real_t distance = position_helper<dims>::distance(data.agent_data.positions.data() + i * dims,
															  data.agent_data.positions.data() + j * dims);

			if (distance <= data.geometries.radius[i] * 2.1)
			{
				data.states.neighbors[i].push_back(j);
				data.rest_lengths[i].push_back(data.geometries.radius[i] * 2);
			}
			else if (distance <= data.geometries.radius[i] * 2.1 * std::sqrt(3))
			{
				data.states.neighbors[i].push_back(j);
				data.rest_lengths[i].push_back(data.geometries.radius[i] * 2 * std::sqrt(3));
			}
		}
	}
}

void kelvin_voigt_solver::update_cell_forces(environment& e)
{
	auto& data = get_cell_data(e);
	update_cell_forces_internal<2>(data, e.cell_definitions_count, e);
}

void kelvin_voigt_model::update_cell_forces(environment& e) { kelvin_voigt_solver::update_cell_forces(e); }

void kelvin_voigt_solver::update_cell_neighbors(environment& e)
{
	auto& data = get_cell_data(e);
	update_cell_neighbors_internal<2>(data);
}

void kelvin_voigt_model::update_cell_neighbors(environment& e) { kelvin_voigt_solver::update_cell_neighbors(e); }
