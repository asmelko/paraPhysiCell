#pragma once

#include <BioFVM/solver.h>

#include "environment.h"
#include "original/modules/pathology.h"
#include "original/modules/settings.h"
#include "simulator_measurer.h"
#include "solver/solver.h"

namespace physicell {

class simulator
{
	biofvm::solver diffusion_solver_;
	physicell::solver mechanics_solver_;

	biofvm::index_t mechanics_step_interval_;
	biofvm::index_t phenotype_step_interval_;
	biofvm::index_t full_save_interval_;
	biofvm::index_t svg_save_interval_;

	void save(environment& e, simulator_durations& durations, PhysiCell_Settings& settings,
			  const cell_coloring_funct_t& cell_coloring_function,
			  const substrate_coloring_funct_t& substrate_coloring_function, biofvm::index_t simulation_step);
	void save_full(environment& e, const PhysiCell_Settings& settings, biofvm::index_t simulation_step);
	void save_svg(environment& e, const PhysiCell_Settings& settings,
				  const cell_coloring_funct_t& cell_coloring_function,
				  const substrate_coloring_funct_t& substrate_coloring_function, biofvm::index_t simulation_step);

	void sync_data_host(environment& e, simulator_durations& durations, biofvm::solvers::data_residency& residency);
	void sync_data_device(environment& e, simulator_durations& durations, biofvm::solvers::data_residency& residency);

public:
	void initialize(environment& e, PhysiCell_Settings& settings);

	void simulate_diffusion(environment& e, simulator_durations& durations, bool& recompute_secretion_and_uptake,
							biofvm::solvers::data_residency& residency);
	void simulate_mechanics(environment& e, simulator_durations& durations, bool& recompute_secretion_and_uptake,
							biofvm::solvers::data_residency& residency);
	void simulate_phenotype(environment& e, simulator_durations& durations, bool& recompute_secretion_and_uptake,
							biofvm::solvers::data_residency& residency);

	void run(environment& e, PhysiCell_Settings& settings, cell_coloring_funct_t cell_coloring_function,
			 substrate_coloring_funct_t substrate_coloring_function = nullptr);
};

} // namespace physicell
