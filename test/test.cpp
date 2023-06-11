#include <BioFVM/microenvironment.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "environment.h"

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

environment default_env(microenvironment& m)
{
	index_t cell_defs_count = 3;

	environment e(m);
	e.cell_definitions_count = cell_defs_count;

	m.agents = std::make_unique<cell_container>(e);

	return e;
}

TEST(cell_container, add_and_remove)
{
	cartesian_mesh mesh(1, { 0, 0, 0 }, { 80, 0, 0 }, { 20, 0, 0 });

	auto m = default_microenv(mesh);
	auto e = default_env(m);

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
