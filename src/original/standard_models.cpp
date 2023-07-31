#include "standard_models.h"

#include "../random.h"
#include "../solver/host/solver_helper.h"
#include "constants.h"
#include "rules.h"
#include "signal_behavior.h"

using namespace biofvm;

namespace physicell {

void standard_cell_transformations(cell& cell, environment& e)
{
	if (cell.phenotype.death.dead() == true)
	{
		return;
	}

	for (index_t i = 0; i < e.cell_definitions_count; i++)
	{
		const auto probability = cell.phenotype.cell_transformations.transformation_rates()[i] * e.phenotype_time_step;
		if (random::instance().uniform() <= probability)
		{
			// std::cout << "Transforming from " << pCell->type_name << " to " << cell_definitions_by_index[i]->name <<
			// std::endl;
			cell.convert(e.cell_definitions[i], i);
			break;
		}
	}
}

void advance_bundled_phenotype_functions(environment& e)
{
	for (auto& cell : e.cast_container<cell_container>().agents())
	{
		// New March 2022
		// perform transformations
		standard_cell_transformations(*cell, e.phenotype_time_step);

		// New March 2023 in Version 1.12.0
		// call the rules-based code to update the phenotype
		if (e.settings.rules_enabled)
		{
			apply_ruleset(cell.get(), e);
		}
		if (get_single_signal(cell.get(), "necrotic") > 0.5)
		{
			real_t rupture = cell->phenotype.volume.rupture_volume();
			real_t volume = cell->phenotype.volume.total();
			if (volume > rupture)
			{
				std::cout << cell->phenotype.volume.total() << " vs " << cell->phenotype.volume.rupture_volume()
						  << " dead: " << get_single_signal(cell.get(), "dead") << std::endl;
				std::cout << cell->phenotype.cycle.current_phase_index() << " "
						  << cell->phenotype.cycle.pCycle_Model->name << std::endl;
			}
		}

		// call the custom code to update the phenotype
		if (cell->functions.update_phenotype)
		{
			cell->functions.update_phenotype(*cell, e.phenotype_time_step);
		}

		// update volume
		if (cell->functions.volume_update_function)
		{
			cell->functions.volume_update_function(*cell, e.phenotype_time_step);
		}

		// update geometry
		cell->phenotype.geometry.update();

		// check for new death events
		if (cell->phenotype.death.check_for_death(e.phenotype_time_step) == true)
		{
			// if so, change the cycle model to the current death model
			cell->phenotype.cycle.sync_to_cycle_model(cell->phenotype.death.current_model());

			// also, turn off motility.
			cell->phenotype.motility.is_motile() = false;

			// turn off secretion, and reduce uptake by a factor of 10
			cell->phenotype.secretion.set_all_secretion_to_zero();
			cell->phenotype.secretion.scale_all_uptake_by_factor(0.10);

			// make sure to run the death entry function
			if (cell->phenotype.cycle.current_phase().entry_function)
			{
				cell->phenotype.cycle.current_phase().entry_function(*cell, e.phenotype_time_step);
			}
		}

		// advance cycle model (for both cell cycle and death cycle models)
		cell->phenotype.cycle.advance_cycle(*cell, e.phenotype_time_step);
	}
}

void chemotaxis_function(cell& cell)
{
	// bias direction is gradient for the indicated substrate
	auto g = cell.nearest_gradient(cell.phenotype.motility.chemotaxis_index());

	for (index_t i = 0; i < cell.e().mechanics_mesh.dims; i++)
	{
		// move up or down gradient based on this direction
		cell.phenotype.motility.migration_bias_direction()[i] = g[i] * cell.phenotype.motility.chemotaxis_direction();
	}

	// normalize
	normalize(cell.phenotype.motility.migration_bias_direction(), cell.e().mechanics_mesh.dims);

	return;
}

template <bool do_normalize>
void advanced_chemotaxis_function(cell& cell)
{
	for (index_t i = 0; i < cell.e().mechanics_mesh.dims; i++)
	{
		cell.phenotype.motility.migration_bias_direction()[i] = 0;
	}

	// weighted combination of the gradients
	for (index_t i = 0; i < cell.e().m.substrates_count; i++)
	{
		// get and normalize ith gradient
		auto g = cell.nearest_gradient(cell.phenotype.motility.chemotaxis_index());

		if constexpr (do_normalize)
			normalize(g.data(), cell.e().mechanics_mesh.dims);

		for (index_t d = 0; d < cell.e().mechanics_mesh.dims; d++)
		{
			cell.phenotype.motility.migration_bias_direction()[d] +=
				g[d] * cell.phenotype.motility.chemotactic_sensitivities()[i];
		}
	}

	// normalize
	normalize(cell.phenotype.motility.migration_bias_direction(), cell.e().mechanics_mesh.dims);
}

void standard_volume_update_function(cell& cell)
{
	cell.phenotype.volume.fluid() += cell.e().phenotype_time_step * cell.phenotype.volume.fluid_change_rate()
									 * (cell.phenotype.volume.target_fluid_fraction() * cell.phenotype.volume.total()
										- cell.phenotype.volume.fluid());

	// if the fluid volume is negative, set to zero
	if (cell.phenotype.volume.fluid() < 0.0)
	{
		cell.phenotype.volume.fluid() = 0.0;
	}

	cell.phenotype.volume.nuclear_fluid() =
		(cell.phenotype.volume.nuclear() / (cell.phenotype.volume.total() + 1e-16)) * (cell.phenotype.volume.fluid());
	cell.phenotype.volume.cytoplasmic_fluid() = cell.phenotype.volume.fluid() - cell.phenotype.volume.nuclear_fluid();

	cell.phenotype.volume.nuclear_solid() +=
		cell.e().phenotype_time_step * cell.phenotype.volume.nuclear_biomass_change_rate()
		* (cell.phenotype.volume.target_solid_nuclear() - cell.phenotype.volume.nuclear_solid());
	if (cell.phenotype.volume.nuclear_solid() < 0.0)
	{
		cell.phenotype.volume.nuclear_solid() = 0.0;
	}

	cell.phenotype.volume.target_solid_cytoplasmic() = cell.phenotype.volume.target_cytoplasmic_to_nuclear_ratio()
													   * // cell.phenotype.volume.cytoplasmic_to_nuclear_fraction() *
													   cell.phenotype.volume.target_solid_nuclear();

	cell.phenotype.volume.cytoplasmic_solid() +=
		cell.e().phenotype_time_step * cell.phenotype.volume.cytoplasmic_biomass_change_rate()
		* (cell.phenotype.volume.target_solid_cytoplasmic() - cell.phenotype.volume.cytoplasmic_solid());
	if (cell.phenotype.volume.cytoplasmic_solid() < 0.0)
	{
		cell.phenotype.volume.cytoplasmic_solid() = 0.0;
	}

	cell.phenotype.volume.solid() = cell.phenotype.volume.nuclear_solid() + cell.phenotype.volume.cytoplasmic_solid();

	cell.phenotype.volume.nuclear() = cell.phenotype.volume.nuclear_solid() + cell.phenotype.volume.nuclear_fluid();
	cell.phenotype.volume.cytoplasmic() =
		cell.phenotype.volume.cytoplasmic_solid() + cell.phenotype.volume.cytoplasmic_fluid();

	cell.phenotype.volume.calcified_fraction() += cell.e().phenotype_time_step
												  * cell.phenotype.volume.calcification_rate()
												  * (1 - cell.phenotype.volume.calcified_fraction());

	cell.phenotype.volume.total() = cell.phenotype.volume.cytoplasmic() + cell.phenotype.volume.nuclear();

	cell.phenotype.volume.fluid_fraction() = cell.phenotype.volume.fluid() / (1e-16 + cell.phenotype.volume.total());

	cell.phenotype.geometry.update();

	return;
}

void standard_elastic_contract_function(cell& lhs, cell& rhs)
{
	const real_t adhesion = sqrt(lhs.phenotype.mechanics.attachment_elastic_constant()
								 * rhs.phenotype.mechanics.attachment_elastic_constant()
								 * lhs.phenotype.mechanics.cell_adhesion_affinities()[rhs.cell_definition_index()]
								 * rhs.phenotype.mechanics.cell_adhesion_affinities()[lhs.cell_definition_index()]);

	for (index_t d = 0; d < lhs.e().mechanics_mesh.dims; d++)
	{
		lhs.velocity()[d] = adhesion * (lhs.position()[d] - rhs.position()[d]);
	}
}

void evaluate_interactions(environment& e)
{
	auto& cells = e.cast_container<cell_container>();

	for (auto& cell : cells.agents())
	{
		if (cell->functions.contact_function)
		{
			for (auto& other : cell->state.attached_cells())
			{
				auto& other_cell = *cells.agents()[other];

				cell->functions.contact_function(*cell, other_cell);
			}
		}
	}
}

void update_cell_and_death_parameters_O2_based(cell& cell, real_t)
{
	// supported cycle models:
	// advanced_Ki67_cycle_model= 0;
	// basic_Ki67_cycle_model=1
	// live_cells_cycle_model = 5;

	if (cell.phenotype.death.dead())
	{
		return;
	}

	// set up shortcuts to find the Q and K(1) phases (assuming Ki67 basic or advanced model)
	static bool indices_initiated = false;
	static int start_phase_index; // Q_phase_index;
	static int end_phase_index;	  // K_phase_index;
	static int necrosis_index;

	static int oxygen_substrate_index = cell.e().m.find_substrate_index("oxygen");

	if (indices_initiated == false)
	{
		// Ki67 models

		if (cell.phenotype.cycle.model().code == constants::advanced_Ki67_cycle_model
			|| cell.phenotype.cycle.model().code == constants::basic_Ki67_cycle_model)
		{
			start_phase_index = cell.phenotype.cycle.model().find_phase_index(constants::Ki67_negative);
			necrosis_index = cell.phenotype.death.find_death_model_index(constants::necrosis_death_model);

			if (cell.phenotype.cycle.model().code == constants::basic_Ki67_cycle_model)
			{
				end_phase_index = cell.phenotype.cycle.model().find_phase_index(constants::Ki67_positive);
				indices_initiated = true;
			}
			if (cell.phenotype.cycle.model().code == constants::advanced_Ki67_cycle_model)
			{
				end_phase_index = cell.phenotype.cycle.model().find_phase_index(constants::Ki67_positive_premitotic);
				indices_initiated = true;
			}
		}

		// live model

		if (cell.phenotype.cycle.model().code == constants::live_cells_cycle_model)
		{
			start_phase_index = cell.phenotype.cycle.model().find_phase_index(constants::live);
			necrosis_index = cell.phenotype.death.find_death_model_index(constants::necrosis_death_model);
			end_phase_index = cell.phenotype.cycle.model().find_phase_index(constants::live);
			indices_initiated = true;
		}

		// cytometry models

		if (cell.phenotype.cycle.model().code == constants::flow_cytometry_cycle_model
			|| cell.phenotype.cycle.model().code == constants::flow_cytometry_separated_cycle_model)
		{
			start_phase_index = cell.phenotype.cycle.model().find_phase_index(constants::G0G1_phase);
			necrosis_index = cell.phenotype.death.find_death_model_index(constants::necrosis_death_model);
			end_phase_index = cell.phenotype.cycle.model().find_phase_index(constants::S_phase);
			indices_initiated = true;
		}

		if (cell.phenotype.cycle.model().code == constants::cycling_quiescent_model)
		{
			start_phase_index = cell.phenotype.cycle.model().find_phase_index(constants::quiescent);
			necrosis_index = cell.phenotype.death.find_death_model_index(constants::necrosis_death_model);
			end_phase_index = cell.phenotype.cycle.model().find_phase_index(constants::cycling);
			indices_initiated = true;
		}
	}

	// don't continue if we never "figured out" the current cycle model.
	if (indices_initiated == false)
	{
		return;
	}

	// sample the microenvironment to get the pO2 value

	double pO2 = (cell.nearest_density_vector())[oxygen_substrate_index]; // constants::oxygen_index];

	// this multiplier is for linear interpolation of the oxygen value
	double multiplier = 1.0;
	if (pO2 < cell.parameters.o2_proliferation_saturation)
	{
		multiplier = (pO2 - cell.parameters.o2_proliferation_threshold)
					 / (cell.parameters.o2_proliferation_saturation - cell.parameters.o2_proliferation_threshold);
	}
	if (pO2 < cell.parameters.o2_proliferation_threshold)
	{
		multiplier = 0.0;
	}

	// now, update the appropriate cycle transition rate

	cell.phenotype.cycle.data.transition_rate(start_phase_index, end_phase_index) =
		multiplier
		* cell.parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index, end_phase_index);

	// Update necrosis rate
	multiplier = 0.0;
	if (pO2 < cell.parameters.o2_necrosis_threshold)
	{
		multiplier = (cell.parameters.o2_necrosis_threshold - pO2)
					 / (cell.parameters.o2_necrosis_threshold - cell.parameters.o2_necrosis_max);
	}
	if (pO2 < cell.parameters.o2_necrosis_max)
	{
		multiplier = 1.0;
	}

	// now, update the necrosis rate

	cell.phenotype.death.rates[necrosis_index] = multiplier * cell.parameters.max_necrosis_rate;

	// check for deterministic necrosis

	if (cell.parameters.necrosis_type == constants::deterministic_necrosis && multiplier > 1e-16)
	{
		cell.phenotype.death.rates[necrosis_index] = 9e99;
	}

	return;
}

} // namespace physicell
