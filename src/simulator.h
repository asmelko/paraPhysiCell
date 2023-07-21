#pragma once

#include <BioFVM/solver.h>

#include "environment.h"
#include "solver/solver.h"

namespace physicell {

class simulator
{
	biofvm::solver diffusion_solver_;
	physicell::solver mechanics_solver_;

    biofvm::index_t simulation_step_;

    biofvm::index_t mechanics_step_interval_;
    biofvm::index_t phenotype_step_interval_;

	bool recompute_secretion_and_uptake_;

public:
	void initialize(environment& e);

    // Called once each diffusion time ticks 
	void simulate_diffusion_and_mechanics(environment& e);
};

} // namespace physicell
