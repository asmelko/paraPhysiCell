#include <BioFVM/microenvironment.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "cell.h"
#include "environment.h"
#include "solver/host/mechanics_solver.h"
#include "solver/host/position_solver.h"

using namespace biofvm;
using namespace physicell;

microenvironment default_microenv(cartesian_mesh mesh)
{
	real_t diffusion_time_step = 5;
	index_t substrates_count = 2;

	auto diff_coefs = std::make_unique<real_t[]>(2);
	diff_coefs[0] = 4;
	diff_coefs[1] = 2;
	auto decay_rates = std::make_unique<real_t[]>(2);
	decay_rates[0] = 5;
	decay_rates[1] = 3;

	auto initial_conds = std::make_unique<real_t[]>(2);
	initial_conds[0] = 1;
	initial_conds[1] = 1;

	microenvironment m(mesh, substrates_count, diffusion_time_step, initial_conds.get());
	m.diffustion_coefficients = std::move(diff_coefs);
	m.decay_rates = std::move(decay_rates);

	return m;
}

environment default_env(microenvironment& m, cartesian_mesh mesh)
{
	index_t cell_defs_count = 3;

	environment e(m, mesh);
	e.cell_definitions_count = cell_defs_count;

	m.agents = std::make_unique<cell_container>(e);

	return e;
}

TEST(cell_container, add_and_remove)
{
	cartesian_mesh mesh(1, { 0, 0, 0 }, { 80, 0, 0 }, { 20, 0, 0 });

	auto m = default_microenv(mesh);
	auto e = default_env(m, mesh);

	cell_container& cont = e.cast_container<cell_container>();

	auto c1 = cont.create();
	c1->volume() = 1;
	c1->velocity()[0] = 10;

	auto c2 = cont.create();
	c2->volume() = 2;
	c2->velocity()[0] = 20;

	auto c3 = cont.create();
	c3->volume() = 3;
	c3->velocity()[0] = 30;

	auto c4 = cont.create();
	c4->volume() = 4;
	c4->velocity()[0] = 40;

	// remove one to last
	{
		cont.remove_agent(c3);

		EXPECT_EQ(cont.data().agent_data.volumes[0], 1);
		EXPECT_EQ(cont.data().agent_data.volumes[1], 2);
		EXPECT_EQ(cont.data().agent_data.volumes[2], 4);

		EXPECT_EQ(cont.data().velocities[0], 10);
		EXPECT_EQ(cont.data().velocities[1], 20);
		EXPECT_EQ(cont.data().velocities[2], 40);
	}

	// remove first
	{
		cont.remove_agent(c1);

		EXPECT_EQ(cont.data().agent_data.volumes[0], 4);
		EXPECT_EQ(cont.data().agent_data.volumes[1], 2);

		EXPECT_EQ(cont.data().velocities[0], 40);
		EXPECT_EQ(cont.data().velocities[1], 20);
	}

	// add one
	{
		auto c5 = cont.create();
		c5->volume() = 5;
		c5->velocity()[0] = 50;

		EXPECT_EQ(cont.data().agent_data.volumes[0], 4);
		EXPECT_EQ(cont.data().agent_data.volumes[1], 2);
		EXPECT_EQ(cont.data().agent_data.volumes[2], 5);

		EXPECT_EQ(cont.data().velocities[0], 40);
		EXPECT_EQ(cont.data().velocities[1], 20);
		EXPECT_EQ(cont.data().velocities[2], 50);
	}
}

void update_cell_for_velocity(cell* c, index_t offset, index_t dims, index_t cell_def_index)
{
	c->cell_definition_index() = cell_def_index;
	c->is_movable() = true;

	for (index_t i = 0; i < dims; ++i)
	{
		c->position()[i] = 1 + offset + i;
	}

	c->phenotype.geometry.radius() = 4 + offset;
	c->phenotype.mechanics.cell_cell_adhesion_strength() = 5 + offset;
	c->phenotype.mechanics.cell_cell_repulsion_strength() = 6 + offset;
	c->phenotype.mechanics.relative_maximum_adhesion_distance() = 7 + offset;
	c->phenotype.mechanics.cell_adhesion_affinities()[0] = 8 + offset;
	c->phenotype.mechanics.cell_adhesion_affinities()[1] = 9 + offset;
	c->phenotype.mechanics.cell_adhesion_affinities()[2] = 10 + offset;
}

class host_velocity_solver : public testing::TestWithParam<index_t>
{};

INSTANTIATE_TEST_SUITE_P(dims, host_velocity_solver, testing::Values(1, 2, 3));

TEST_P(host_velocity_solver, simple)
{
	auto dims = GetParam();

	cartesian_mesh mesh(dims, { 0, 0, 0 }, { 80, 80, 80 }, { 20, 20, 20 });

	auto m = default_microenv(mesh);
	auto e = default_env(m, mesh);

	cell_container& cont = e.cast_container<cell_container>();

	auto c1 = cont.create();
	update_cell_for_velocity(c1, 0, dims, 0);

	auto c2 = cont.create();
	update_cell_for_velocity(c2, 10, dims, 1);

	auto c3 = cont.create();
	for (index_t i = 0; i < dims; ++i)
		c3->position()[i] = 70;

	mechanics_solver::update_mechanics_mesh(e);
	position_solver::update_cell_velocities_and_neighbors(e);

	real_t expected;

	if (dims == 3)
		expected = 55.61362893236258;
	else if (dims == 2)
		expected = 69.556601146394243;
	else
		expected = 100.15967659502964;

	for (index_t i = 0; i < dims; ++i)
	{
		EXPECT_FLOAT_EQ(c1->velocity()[i], expected);
		EXPECT_FLOAT_EQ(c2->velocity()[i], -expected);
	}

	for (index_t i = 0; i < dims; ++i)
	{
		EXPECT_FLOAT_EQ(c3->velocity()[i], 0);
	}
}

std::vector<real_t> compute_expected_velocities(const std::vector<std::unique_ptr<cell>>& cells, index_t dims)
{
	std::vector<real_t> expected_velocities(cells.size() * dims, 0);

	for (std::size_t i = 0; i < cells.size(); ++i)
	{
		auto& c1 = cells[i];

		for (std::size_t j = 0; j < cells.size(); j++)
		{
			if (i == j)
				continue;

			auto& c2 = cells[j];

			std::vector<real_t> diff(dims, 0);
			real_t dist = 0;

			for (index_t k = 0; k < dims; ++k)
			{
				diff[k] = c1->position()[k] - c2->position()[k];
				dist += diff[k] * diff[k];
			}

			dist = std::sqrt(dist);

			real_t repulsion;
			{
				const real_t repulsive_distance = c1->phenotype.geometry.radius() + c2->phenotype.geometry.radius();

				repulsion = 1 - dist / repulsive_distance;

				repulsion *= repulsion;

				repulsion *= sqrt(c1->phenotype.mechanics.cell_cell_repulsion_strength()
								  * c2->phenotype.mechanics.cell_cell_repulsion_strength());
			}

			// compute adhesion
			real_t adhesion;
			{
				const real_t adhesion_distance =
					c1->phenotype.mechanics.relative_maximum_adhesion_distance() * c1->phenotype.geometry.radius()
					+ c2->phenotype.mechanics.relative_maximum_adhesion_distance() * c2->phenotype.geometry.radius();

				adhesion = 1 - dist / adhesion_distance;

				adhesion *= adhesion;

				const index_t lhs = c1->cell_definition_index();
				const index_t rhs = c2->cell_definition_index();

				const real_t lhs_cell_adhesion_affinity = c1->phenotype.mechanics.cell_adhesion_affinities()[rhs];
				const real_t rhs_cell_adhesion_affinity = c2->phenotype.mechanics.cell_adhesion_affinities()[lhs];

				adhesion *= sqrt(c1->phenotype.mechanics.cell_cell_adhesion_strength()
								 * c2->phenotype.mechanics.cell_cell_adhesion_strength() * lhs_cell_adhesion_affinity
								 * rhs_cell_adhesion_affinity);
			}

			real_t force = (repulsion - adhesion) / dist;

			for (index_t k = 0; k < dims; k++)
			{
				expected_velocities[i * dims + k] += force * diff[k];
			}
		}
	}

	return expected_velocities;
}

TEST_P(host_velocity_solver, complex)
{
	auto dims = GetParam();

	cartesian_mesh mesh(dims, { 0, 0, 0 }, { 80, 80, 80 }, { 20, 20, 20 });

	auto m = default_microenv(mesh);
	auto e = default_env(m, mesh);

	cell_container& cont = e.cast_container<cell_container>();

	auto c1 = cont.create();
	update_cell_for_velocity(c1, 0, dims, 0);
	c1->phenotype.geometry.radius() = dims == 1 ? 2 : 4;
	c1->phenotype.mechanics.relative_maximum_adhesion_distance() = 1;
	for (index_t i = 0; i < dims; ++i)
		c1->position()[i] = 1;

	auto c2 = cont.create();
	update_cell_for_velocity(c2, 10, dims, 1);
	c2->phenotype.geometry.radius() = dims == 1 ? 2 : 4;
	c2->phenotype.mechanics.relative_maximum_adhesion_distance() = 1;
	for (index_t i = 0; i < dims; ++i)
		c2->position()[i] = 4;

	auto c3 = cont.create();
	c3->phenotype.geometry.radius() = dims == 1 ? 2 : 4;
	c3->phenotype.mechanics.relative_maximum_adhesion_distance() = 1;
	for (index_t i = 0; i < dims; ++i)
		c3->position()[i] = 7;

	mechanics_solver::update_mechanics_mesh(e);
	position_solver::update_cell_velocities_and_neighbors(e);

	EXPECT_EQ(cont.data().neighbors[0].size(), 1);
	EXPECT_EQ(cont.data().neighbors[1].size(), 2);
	EXPECT_EQ(cont.data().neighbors[2].size(), 1);

	auto expected_velocities = compute_expected_velocities(cont.agents(), dims);

	for (index_t i = 0; i < dims; i++)
		EXPECT_FLOAT_EQ(c1->velocity()[i], expected_velocities[i]);

	for (index_t i = 0; i < dims; i++)
		EXPECT_FLOAT_EQ(c2->velocity()[i], expected_velocities[dims + i]);

	for (index_t i = 0; i < dims; i++)
		EXPECT_FLOAT_EQ(c3->velocity()[i], expected_velocities[2 * dims + i]);
}
