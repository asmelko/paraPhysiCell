#pragma once

#include <chrono>

#include <BioFVM/solver.h>

#include "environment.h"
#include "original/modules/pathology.h"
#include "original/modules/settings.h"
#include "solver/solver.h"

namespace physicell {

struct simulator_durations
{
	std::chrono::high_resolution_clock::time_point start, end;

	std::size_t diffusion = 0;
	std::size_t secretion = 0;
	std::size_t gradient = 0;
	std::size_t custom_rules = 0;
	std::size_t custom_interactions = 0;
	std::size_t forces = 0;
	std::size_t motility = 0;
	std::size_t membrane = 0;
	std::size_t spring = 0;
	std::size_t position = 0;
	std::size_t cell_cell = 0;
	std::size_t container_mech = 0;
	std::size_t mesh_mech = 0;
	std::size_t neighbors_mech = 0;
	std::size_t advance_phe = 0;
	std::size_t container_phe = 0;
	std::size_t mesh_phe = 0;
	std::size_t neighbors_phe = 0;
	std::size_t full_save = 0;
	std::size_t svg_save = 0;

	void print_durations();
};

class simulator
{
	simulator_durations durations_;

	biofvm::solver diffusion_solver_;
	physicell::solver mechanics_solver_;

	biofvm::index_t simulation_step_;

	biofvm::index_t mechanics_step_interval_;
	biofvm::index_t phenotype_step_interval_;
	biofvm::index_t full_save_interval_;
	biofvm::index_t svg_save_interval_;

	bool recompute_secretion_and_uptake_;

	void save_full(environment& e, PhysiCell_Settings& settings);
	void save_svg(environment& e, PhysiCell_Settings& settings, const cell_coloring_funct_t& cell_coloring_function);

public:
	void initialize(environment& e, PhysiCell_Settings& settings);

	// Called once each diffusion time ticks
	void simulate_diffusion_and_mechanics(environment& e);

	void run(environment& e, PhysiCell_Settings& settings, cell_coloring_funct_t cell_coloring_function);
};

} // namespace physicell
