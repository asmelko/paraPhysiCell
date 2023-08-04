#pragma once

#include <cmath>

#include <BioFVM/mesh.h>

#include "../../random.h"

namespace physicell {

constexpr biofvm::real_t zero_threshold = 1e-16;

constexpr void update_membrane_velocity(biofvm::real_t position, biofvm::real_t bounding_box, biofvm::real_t sign,
										biofvm::real_t radius, biofvm::real_t repulsion_strength,
										biofvm::real_t& velocity)
{
	biofvm::real_t distance = std::abs(bounding_box - position);

	distance = std::max<biofvm::real_t>(distance, 0.00001);

	biofvm::real_t repulsion = 1 - distance / radius;
	repulsion = repulsion < 0 ? 0 : repulsion;

	repulsion *= repulsion * repulsion_strength * sign;

	velocity += repulsion * distance;
}

constexpr void normalize(biofvm::real_t* __restrict__ vector, biofvm::index_t dims)
{
	biofvm::real_t length = 0;

	for (biofvm::index_t d = 0; d < dims; d++)
	{
		length += vector[d] * vector[d];
	}

	length = std::sqrt(length);

	for (biofvm::index_t d = 0; d < dims; d++)
	{
		vector[d] = length > zero_threshold ? vector[d] / length : 0;
	}
}

template <biofvm::index_t dims>
struct position_helper
{};

template <>
struct position_helper<1>
{
	static constexpr biofvm::real_t difference_and_distance(const biofvm::real_t* __restrict__ lhs,
															const biofvm::real_t* __restrict__ rhs,
															biofvm::real_t* __restrict__ difference)
	{
		difference[0] = lhs[0] - rhs[0];

		return std::abs(difference[0]);
	}

	static constexpr void update_velocities(biofvm::real_t* __restrict__ lhs, biofvm::real_t* __restrict__ rhs,
											const biofvm::real_t* __restrict__ difference, const biofvm::real_t force)
	{
		lhs[0] += force * difference[0];
		rhs[0] -= force * difference[0];
	}

	static void random_walk(bool, biofvm::real_t* __restrict__ walk)
	{
		biofvm::real_t rand = random::instance().uniform();
		walk[0] = rand < 0.5 ? -1 : 1;
	}

	static constexpr void update_motility_vector(biofvm::real_t* __restrict__ motility_vector,
												 const biofvm::real_t* __restrict__ walk,
												 const biofvm::real_t* __restrict__ migration_bias_direction,
												 const biofvm::real_t migration_bias)
	{
		motility_vector[0] = (1 - migration_bias) * walk[0] + migration_bias * migration_bias_direction[0];
	}

	static constexpr void normalize_and_scale(biofvm::real_t* __restrict__ vector, biofvm::real_t scale)
	{
		biofvm::real_t length = std::abs(vector[0]);

		vector[0] = length > zero_threshold ? vector[0] * scale / length : 0;
	}

	static constexpr void update_membrane_velocities(biofvm::real_t* __restrict__ velocity,
													 const biofvm::real_t* __restrict__ position,
													 const biofvm::cartesian_mesh& mesh, const biofvm::real_t radius,
													 const biofvm::real_t repulsion_strength)
	{
		update_membrane_velocity(position[0], mesh.bounding_box_mins[0], 1, radius, repulsion_strength, velocity[0]);
		update_membrane_velocity(position[0], mesh.bounding_box_maxs[1], -1, radius, repulsion_strength, velocity[0]);
	}

	static constexpr void add(biofvm::real_t* __restrict__ lhs, const biofvm::real_t* __restrict__ rhs)
	{
		lhs[0] += rhs[0];
	}

	static constexpr void subtract(biofvm::real_t* __restrict__ dst, const biofvm::real_t* __restrict__ lhs,
								   const biofvm::real_t* __restrict__ rhs)
	{
		dst[0] = lhs[0] - rhs[0];
	}
};

template <>
struct position_helper<2>
{
	static constexpr biofvm::real_t difference_and_distance(const biofvm::real_t* __restrict__ lhs,
															const biofvm::real_t* __restrict__ rhs,
															biofvm::real_t* __restrict__ difference)
	{
		difference[0] = lhs[0] - rhs[0];
		difference[1] = lhs[1] - rhs[1];

		return std::sqrt(difference[0] * difference[0] + difference[1] * difference[1]);
	}

	static constexpr void update_velocities(biofvm::real_t* __restrict__ lhs, biofvm::real_t* __restrict__ rhs,
											const biofvm::real_t* __restrict__ difference, const biofvm::real_t force)
	{
		lhs[0] += force * difference[0];
		lhs[1] += force * difference[1];

		rhs[0] -= force * difference[0];
		rhs[1] -= force * difference[1];
	}

	static void random_walk(bool, biofvm::real_t* __restrict__ walk)
	{
		biofvm::real_t theta = random::instance().uniform(0, 2 * M_PI);
		walk[0] = std::cos(theta);
		walk[1] = std::sin(theta);
	}

	static constexpr void update_motility_vector(biofvm::real_t* __restrict__ motility_vector,
												 const biofvm::real_t* __restrict__ walk,
												 const biofvm::real_t* __restrict__ migration_bias_direction,
												 const biofvm::real_t migration_bias)
	{
		motility_vector[0] = (1 - migration_bias) * walk[0] + migration_bias * migration_bias_direction[0];
		motility_vector[1] = (1 - migration_bias) * walk[1] + migration_bias * migration_bias_direction[1];
	}

	static constexpr void normalize_and_scale(biofvm::real_t* __restrict__ vector, biofvm::real_t scale)
	{
		biofvm::real_t length = std::sqrt(vector[0] * vector[0] + vector[1] * vector[1]);

		vector[0] = length > zero_threshold ? vector[0] * scale / length : 0;
		vector[1] = length > zero_threshold ? vector[1] * scale / length : 0;
	}

	static constexpr void update_membrane_velocities(biofvm::real_t* __restrict__ velocity,
													 const biofvm::real_t* __restrict__ position,
													 const biofvm::cartesian_mesh& mesh, const biofvm::real_t radius,
													 const biofvm::real_t repulsion_strength)
	{
		update_membrane_velocity(position[0], mesh.bounding_box_mins[0], 1, radius, repulsion_strength, velocity[0]);
		update_membrane_velocity(position[0], mesh.bounding_box_maxs[0], -1, radius, repulsion_strength, velocity[0]);
		update_membrane_velocity(position[1], mesh.bounding_box_mins[1], 1, radius, repulsion_strength, velocity[1]);
		update_membrane_velocity(position[1], mesh.bounding_box_maxs[1], -1, radius, repulsion_strength, velocity[1]);
	}

	static constexpr void add(biofvm::real_t* __restrict__ lhs, const biofvm::real_t* __restrict__ rhs)
	{
		lhs[0] += rhs[0];
		lhs[1] += rhs[1];
	}

	static constexpr void subtract(biofvm::real_t* __restrict__ dst, const biofvm::real_t* __restrict__ lhs,
								   const biofvm::real_t* __restrict__ rhs)
	{
		dst[0] = lhs[0] - rhs[0];
		dst[1] = lhs[1] - rhs[1];
	}
};

template <>
struct position_helper<3>
{
	static constexpr biofvm::real_t difference_and_distance(const biofvm::real_t* __restrict__ lhs,
															const biofvm::real_t* __restrict__ rhs,
															biofvm::real_t* __restrict__ difference)
	{
		difference[0] = lhs[0] - rhs[0];
		difference[1] = lhs[1] - rhs[1];
		difference[2] = lhs[2] - rhs[2];

		return std::sqrt(difference[0] * difference[0] + difference[1] * difference[1] + difference[2] * difference[2]);
	}

	static constexpr void update_velocities(biofvm::real_t* __restrict__ lhs, biofvm::real_t* __restrict__ rhs,
											const biofvm::real_t* __restrict__ difference, const biofvm::real_t force)
	{
		lhs[0] += force * difference[0];
		lhs[1] += force * difference[1];
		lhs[2] += force * difference[2];

		rhs[0] -= force * difference[0];
		rhs[1] -= force * difference[1];
		rhs[2] -= force * difference[2];
	}

	static void random_walk(bool restrict_to_2d, biofvm::real_t* __restrict__ walk)
	{
		if (restrict_to_2d)
		{
			position_helper<2>::random_walk(true, walk);
			walk[2] = 0;
		}
		else
		{
			const biofvm::real_t theta = random::instance().uniform(0, 2 * M_PI);
			const biofvm::real_t z = random::instance().uniform(-1, 1);
			const biofvm::real_t r = std::sqrt(1 - z * z);

			walk[0] = std::cos(theta) * r;
			walk[1] = std::sin(theta) * r;
			walk[2] = z;
		}
	}

	static constexpr void update_motility_vector(biofvm::real_t* __restrict__ motility_vector,
												 const biofvm::real_t* __restrict__ walk,
												 const biofvm::real_t* __restrict__ migration_bias_direction,
												 const biofvm::real_t migration_bias)
	{
		motility_vector[0] = (1 - migration_bias) * walk[0] + migration_bias * migration_bias_direction[0];
		motility_vector[1] = (1 - migration_bias) * walk[1] + migration_bias * migration_bias_direction[1];
		motility_vector[2] = (1 - migration_bias) * walk[2] + migration_bias * migration_bias_direction[2];
	}

	static constexpr void normalize_and_scale(biofvm::real_t* __restrict__ vector, biofvm::real_t scale)
	{
		biofvm::real_t length = std::sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);

		vector[0] = length > zero_threshold ? vector[0] * scale / length : 0;
		vector[1] = length > zero_threshold ? vector[1] * scale / length : 0;
		vector[2] = length > zero_threshold ? vector[2] * scale / length : 0;
	}

	static constexpr void update_membrane_velocities(biofvm::real_t* __restrict__ velocity,
													 const biofvm::real_t* __restrict__ position,
													 const biofvm::cartesian_mesh& mesh, const biofvm::real_t radius,
													 const biofvm::real_t repulsion_strength)
	{
		update_membrane_velocity(position[0], mesh.bounding_box_mins[0], 1, radius, repulsion_strength, velocity[0]);
		update_membrane_velocity(position[0], mesh.bounding_box_maxs[0], -1, radius, repulsion_strength, velocity[0]);
		update_membrane_velocity(position[1], mesh.bounding_box_mins[1], 1, radius, repulsion_strength, velocity[1]);
		update_membrane_velocity(position[1], mesh.bounding_box_maxs[1], -1, radius, repulsion_strength, velocity[1]);
		update_membrane_velocity(position[2], mesh.bounding_box_mins[2], 1, radius, repulsion_strength, velocity[2]);
		update_membrane_velocity(position[2], mesh.bounding_box_maxs[2], -1, radius, repulsion_strength, velocity[2]);
	}

	static constexpr void add(biofvm::real_t* __restrict__ lhs, const biofvm::real_t* __restrict__ rhs)
	{
		lhs[0] += rhs[0];
		lhs[1] += rhs[1];
		lhs[2] += rhs[2];
	}

	static constexpr void subtract(biofvm::real_t* __restrict__ dst, const biofvm::real_t* __restrict__ lhs,
								   const biofvm::real_t* __restrict__ rhs)
	{
		dst[0] = lhs[0] - rhs[0];
		dst[1] = lhs[1] - rhs[1];
		dst[2] = lhs[2] - rhs[2];
	}
};
} // namespace physicell
