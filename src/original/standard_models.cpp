#include "standard_models.h"

#include "../random.h"
#include "../solver/host/solver_helper.h"
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

} // namespace physicell
