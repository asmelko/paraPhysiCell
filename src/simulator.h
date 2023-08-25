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
	using full_simulate_func_t =
		std::function<void(environment&, simulator_durations&, bool&, biofvm::solvers::data_residency&)>;
	using simulate_func_t = std::function<void(environment&)>;

	std::vector<std::pair<biofvm::index_t, full_simulate_func_t>> custom_simulations_;

	environment& e;

	biofvm::solver diffusion_solver_;
	physicell::solver mechanics_solver_;

	biofvm::index_t mechanics_step_interval_;
	biofvm::index_t phenotype_step_interval_;
	biofvm::index_t full_save_interval_;
	biofvm::index_t svg_save_interval_;
	biofvm::index_t max_time_;

	void save(simulator_durations& durations, PhysiCell_Settings& settings,
			  const cell_coloring_funct_t& cell_coloring_function,
			  const substrate_coloring_funct_t& substrate_coloring_function, biofvm::index_t simulation_step);
	void save_full(const PhysiCell_Settings& settings, biofvm::index_t simulation_step);
	void save_svg(const PhysiCell_Settings& settings, const cell_coloring_funct_t& cell_coloring_function,
				  const substrate_coloring_funct_t& substrate_coloring_function, biofvm::index_t simulation_step);

	void sync_data_host(simulator_durations& durations, biofvm::solvers::data_residency& residency);
	void sync_data_device(simulator_durations& durations, biofvm::solvers::data_residency& residency);

public:
	simulator(environment& e);

	void initialize(PhysiCell_Settings& settings);

	void simulate_diffusion_device(simulator_durations& durations, bool& recompute_secretion_and_uptake,
								   biofvm::solvers::data_residency& residency);
	void simulate_diffusion(simulator_durations& durations, bool& recompute_secretion_and_uptake,
							biofvm::solvers::data_residency& residency);
	void simulate_mechanics(simulator_durations& durations, bool& recompute_secretion_and_uptake,
							biofvm::solvers::data_residency& residency);
	void simulate_phenotype(simulator_durations& durations, bool& recompute_secretion_and_uptake,
							biofvm::solvers::data_residency& residency);

	void add_simulation_step(biofvm::real_t time_step, simulate_func_t&& simulate_f);

	void run(PhysiCell_Settings& settings, cell_coloring_funct_t cell_coloring_function,
			 substrate_coloring_funct_t substrate_coloring_function = nullptr);
};

} // namespace physicell
