#include "morse_kv_membrane_solver.h"

#include <BioFVM/microenvironment.h>

#include "../../models/morse_kv_membrane_model.h"
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
void update_cell_forces_internal(cell_data& data, environment& e)
{
#pragma omp for
	for (index_t i = 0; i < data.agents_count; i++)
	{
		if (data.is_movable[i] == 0)
			continue;

		if (data.membrane_neighbors[i] == 0)
			continue;

		for (index_t nei_i = 0; nei_i < (index_t)data.states.neighbors[i].size(); nei_i++)
		{
			solve_pair_intra<dims>(i, data.states.neighbors[i][nei_i], data.rest_lengths[i][nei_i],
								   e.mechanics_time_step, data.velocities.data(), data.agent_data.positions.data(),
								   data.spring_constants.data(), data.dissipation_rates.data(),
								   data.previous_velocities.data());
		}
	}
}

bool membrane_updated = false;
std::vector<index_t> membrane_neighbors_count_initial;
std::vector<index_t> membrane_neighbors_count_current;

template <index_t dims>
void find_external_agents(cell_data& data, environment& e)
{
	if (!membrane_updated)
	{
		membrane_updated = true;

		for (index_t i = 0; i < data.agents_count; i++)
		{
			data.membrane_neighbors[i] = 0;

			std::vector<real_t> angles;

			for (index_t j = 0; j < data.agents_count; j++)
			{
				if (i == j || data.cell_residency[i] != data.cell_residency[j])
					continue;

				real_t displacement[dims];

				const real_t distance = position_helper<dims>::difference_and_distance(
					data.agent_data.positions.data() + j * dims, data.agent_data.positions.data() + i * dims,
					displacement);

				if (distance > data.geometries.radius[i] * 2.1 * std::sqrt(2))
					continue;

				real_t cos_alpha = displacement[1] / distance;
				real_t alpha = std::acos(cos_alpha) * (180 / M_PI);
				if (displacement[0] < 0)
					alpha = 360 - alpha;

				angles.push_back(alpha);
			}

			std::cout << "Cell " << i << " angles: ";
			for (const auto& angle : angles)
			{
				std::cout << angle << " ";
			}
			std::cout << std::endl;

			if (angles.size() <= 3)
			{
				data.membrane_neighbors[i] = 1;
				continue;
			}

			std::sort(angles.begin(), angles.end());

			for (size_t k = 1; k < angles.size(); k++)
			{
				if (angles[k] - angles[k - 1] >= 70)
					data.membrane_neighbors[i] = 1;
			}

			if (360 - angles.back() + angles.front() >= 70)
				data.membrane_neighbors[i] = 1;
		}

		for (index_t i = 0; i < data.agents_count; i++)
		{
			if (data.membrane_neighbors[i] == 0)
				continue;

			membrane_neighbors_count_initial.resize(data.cell_residency[i] + 1);
			membrane_neighbors_count_initial[data.cell_residency[i]]++;

			index_t count = 0;

			for (index_t j = 0; j < data.agents_count; j++)
			{
				if (i == j || data.cell_residency[i] != data.cell_residency[j] || data.membrane_neighbors[j] == 0)
					continue;

				const real_t distance = position_helper<dims>::distance(data.agent_data.positions.data() + j * dims,
																		data.agent_data.positions.data() + i * dims);

				if (distance > data.geometries.radius[i] * 2.1)
					continue;

				// if (count >= 2)
				// 	throw std::runtime_error("More than 2 membrane neighbors found");

				data.states.neighbors[i].push_back(j);
				data.rest_lengths[i].push_back(data.geometries.radius[i] * 2);

				count++;
			}
		}

		membrane_neighbors_count_current = membrane_neighbors_count_initial;
	}

	// shrinking membrane
	// for (index_t i = 0; i < data.agents_count; i++)
	// {
	// 	if (data.membrane_neighbors[i] == 0 || data.states.neighbors[i].size() != 2)
	// 		continue;

	// 	index_t nei_i = data.states.neighbors[i][0];
	// 	index_t nei_j = data.states.neighbors[i][1];

	// 	const real_t distance = position_helper<dims>::distance(data.agent_data.positions.data() + nei_i * dims,
	// 															data.agent_data.positions.data() + nei_j * dims);

	// 	if (distance < data.geometries.radius[i] * 2.3)
	// 	{
	// 		data.membrane_neighbors[i] = 0;
	// 		data.states.neighbors[i].clear();
	// 		data.rest_lengths[i].clear();
	// 		auto it = std::find(data.states.neighbors[nei_i].begin(), data.states.neighbors[nei_i].end(), i);
	// 		if (it == data.states.neighbors[nei_i].end())
	// 			throw std::runtime_error("Neighbor not found");
	// 		*it = nei_j;

	// 		it = std::find(data.states.neighbors[nei_j].begin(), data.states.neighbors[nei_j].end(), i);
	// 		if (it == data.states.neighbors[nei_j].end())
	// 			throw std::runtime_error("Neighbor not found");
	// 		*it = nei_i;
	// 	}
	// }

	for (index_t i = 0; i < data.agents_count; i++)
	{
		if (data.membrane_neighbors[i] == 0)
			continue;

		if ((real_t)membrane_neighbors_count_current[data.cell_residency[i]]
				/ (real_t)membrane_neighbors_count_initial[data.cell_residency[i]]
			>= e.membrane_stretching_factors[data.cell_definition_indices[i]])
			continue;

		for (auto j : data.states.neighbors[i])
		{
			const real_t distance = position_helper<dims>::distance(data.agent_data.positions.data() + j * dims,
																	data.agent_data.positions.data() + i * dims);

			if (distance > data.geometries.radius[i] * 2.5)
			{
				for (index_t k = 0; k < data.agents_count; k++)
				{
					if (i == k || j == k || data.cell_residency[i] != data.cell_residency[k]
						|| data.cell_residency[j] != data.cell_residency[k] || data.membrane_neighbors[k] == 1)
						continue;

					const real_t distance_i = position_helper<dims>::distance(
						data.agent_data.positions.data() + k * dims, data.agent_data.positions.data() + i * dims);

					const real_t distance_j = position_helper<dims>::distance(
						data.agent_data.positions.data() + k * dims, data.agent_data.positions.data() + j * dims);

					if (distance_i + distance_j > data.geometries.radius[i] * 4.2)
						continue;

					membrane_neighbors_count_current[data.cell_residency[i]]++;
					data.membrane_neighbors[k] = 1;

					auto it = std::find(data.states.neighbors[i].begin(), data.states.neighbors[i].end(), j);
					if (it == data.states.neighbors[i].end())
						throw std::runtime_error("Neighbor not found");
					*it = k;

					it = std::find(data.states.neighbors[j].begin(), data.states.neighbors[j].end(), i);
					if (it == data.states.neighbors[j].end())
						throw std::runtime_error("Neighbor not found");
					*it = k;

					data.states.neighbors[k].clear();
					data.states.neighbors[k].push_back(i);
					data.states.neighbors[k].push_back(j);
					data.rest_lengths[k].clear();
					data.rest_lengths[k].push_back(data.geometries.radius[k] * 2);
					data.rest_lengths[k].push_back(data.geometries.radius[k] * 2);


					for (index_t k = 0; k < data.agents_count; k++)
					{
						if (data.cell_residency[i] != data.cell_residency[k] || data.membrane_neighbors[k] == 0)
							continue;

						data.spring_constants[k] +=
							(e.membrane_stretched_spring_stepup[data.cell_definition_indices[i]])
							/ (membrane_neighbors_count_initial[data.cell_residency[i]]
							   * (e.membrane_stretching_factors[data.cell_definition_indices[i]] - 1));
					}
				}
			}
		}
	}
}

void morse_kv_membrane_solver::update_cell_forces(environment& e)
{
	auto& data = get_cell_data(e);
	update_cell_forces_internal<2>(data, e);
}

void morse_kv_membrane_position_model::update_cell_forces(environment& e)
{
	morse_kv_membrane_solver::update_cell_forces(e);
}

void morse_kv_membrane_solver::update_cell_neighbors(environment& e)
{
	auto& data = get_cell_data(e);
#pragma omp critical
	find_external_agents<2>(data, e);
#pragma omp barrier
}

void morse_kv_membrane_position_model::update_cell_neighbors(environment& e)
{
	morse_kv_membrane_solver::update_cell_neighbors(e);
}
