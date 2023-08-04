#include "standard_models.h"

#include "../random.h"
#include "../solver/host/solver_helper.h"
#include "constants.h"
#include "rules.h"
#include "signal_behavior.h"

using namespace biofvm;

namespace physicell {

bool PhysiCell_standard_models_initialized = false;
bool PhysiCell_standard_death_models_initialized = false;
bool PhysiCell_standard_cycle_models_initialized = false;

Cycle_Model Ki67_advanced, Ki67_basic, live, apoptosis, necrosis;
Cycle_Model cycling_quiescent;
Death_Parameters apoptosis_parameters, necrosis_parameters;

// new cycle models:

Cycle_Model flow_cytometry_cycle_model, flow_cytometry_separated_cycle_model;

void standard_Ki67_positive_phase_entry_function(cell& cell, real_t)
{
	// the cell wants to double its volume
	cell.phenotype.volume.target_solid_nuclear() *= 2.0;
	cell.phenotype.volume.target_solid_cytoplasmic() *= 2.0;
	return;
}

void standard_Ki67_negative_phase_entry_function(cell&, real_t) { return; }

void standard_live_phase_entry_function(cell& cell, real_t)
{
	// the cell wants to double its volume
	cell.phenotype.volume.target_solid_nuclear() *= 2.0;
	cell.phenotype.volume.target_solid_cytoplasmic() *= 2.0;

	return;
}

void S_phase_entry_function(cell& cell, real_t)
{
	// the cell wants to double its volume
	cell.phenotype.volume.target_solid_nuclear() *= 2.0;
	cell.phenotype.volume.target_solid_cytoplasmic() *= 2.0;

	return;
}

void standard_cycling_entry_function(cell& cell, real_t)
{
	// the cell wants to double its volume
	cell.phenotype.volume.target_solid_nuclear() *= 2.0;
	cell.phenotype.volume.target_solid_cytoplasmic() *= 2.0;

	return;
}

void standard_apoptosis_entry_function(cell& cell, real_t)
{
	// the volume model wants to shrink the cell
	cell.phenotype.volume.target_fluid_fraction() = 0.0;
	cell.phenotype.volume.target_solid_cytoplasmic() = 0.0;
	cell.phenotype.volume.target_solid_nuclear() = 0.0;

	cell.phenotype.volume.target_cytoplasmic_to_nuclear_ratio() = 0.0;

	// change the rate parameters

	cell.phenotype.volume.cytoplasmic_biomass_change_rate() =
		cell.phenotype.death.current_parameters().cytoplasmic_biomass_change_rate;
	cell.phenotype.volume.nuclear_biomass_change_rate() =
		cell.phenotype.death.current_parameters().nuclear_biomass_change_rate;
	cell.phenotype.volume.fluid_change_rate() = cell.phenotype.death.current_parameters().unlysed_fluid_change_rate;

	cell.phenotype.volume.calcification_rate() = cell.phenotype.death.current_parameters().calcification_rate;

	cell.phenotype.volume.relative_rupture_volume() = cell.phenotype.death.current_parameters().relative_rupture_volume;
	cell.phenotype.volume.rupture_volume() =
		cell.phenotype.volume.total() * cell.phenotype.volume.relative_rupture_volume();

	return;
}

void standard_necrosis_entry_function(cell& cell, real_t)
{
	// the volume model wants to degrade the solids, but swell by osmosis
	cell.phenotype.volume.target_fluid_fraction() = 1.0;
	cell.phenotype.volume.target_solid_cytoplasmic() = 0.0;
	cell.phenotype.volume.target_solid_nuclear() = 0.0;

	cell.phenotype.volume.target_cytoplasmic_to_nuclear_ratio() = 0.0;

	// change the rate parameters
	cell.phenotype.volume.cytoplasmic_biomass_change_rate() =
		cell.phenotype.death.current_parameters().cytoplasmic_biomass_change_rate;
	cell.phenotype.volume.nuclear_biomass_change_rate() =
		cell.phenotype.death.current_parameters().nuclear_biomass_change_rate;
	cell.phenotype.volume.fluid_change_rate() = cell.phenotype.death.current_parameters().unlysed_fluid_change_rate;

	cell.phenotype.volume.calcification_rate() = cell.phenotype.death.current_parameters().calcification_rate;

	// set the bursting volume
	cell.phenotype.volume.relative_rupture_volume() = cell.phenotype.death.current_parameters().relative_rupture_volume;
	cell.phenotype.volume.rupture_volume() =
		cell.phenotype.volume.total() * cell.phenotype.volume.relative_rupture_volume();

	return;
}

void standard_lysis_entry_function(cell& cell, real_t)
{
	// the volume model wants to shrink the cell
	cell.phenotype.volume.target_fluid_fraction() = 0.0;
	cell.phenotype.volume.target_solid_cytoplasmic() = 0.0;
	cell.phenotype.volume.target_solid_nuclear() = 0.0;

	// change the rate parameters

	cell.phenotype.volume.cytoplasmic_biomass_change_rate() =
		cell.phenotype.death.current_parameters().cytoplasmic_biomass_change_rate;
	cell.phenotype.volume.nuclear_biomass_change_rate() =
		cell.phenotype.death.current_parameters().nuclear_biomass_change_rate;
	cell.phenotype.volume.fluid_change_rate() = cell.phenotype.death.current_parameters().lysed_fluid_change_rate;

	cell.phenotype.volume.calcification_rate() = cell.phenotype.death.current_parameters().calcification_rate;

	// set the bursting volume
	cell.phenotype.volume.relative_rupture_volume() = 9e99;
	cell.phenotype.volume.rupture_volume() =
		cell.phenotype.volume.total() * cell.phenotype.volume.relative_rupture_volume();

	return;
}

bool standard_necrosis_arrest_function(cell& cell, real_t)
{
	// remain in the non-lysed state / phase if volume has not exceeded the
	// rupture volume
	if (cell.phenotype.volume.total() < cell.phenotype.volume.rupture_volume())
	{
		return true;
	}

	return false;
}

/* create standard models */

void create_ki67_models()
{
	// Ki67_basic:

	Ki67_basic.code = constants::basic_Ki67_cycle_model;
	Ki67_basic.name = "Ki67 (basic)";

	Ki67_basic.data.time_units = "min";

	Ki67_basic.add_phase(constants::Ki67_negative, "Ki67-");
	Ki67_basic.add_phase(constants::Ki67_positive, "Ki67+");

	Ki67_basic.phases[1].division_at_phase_exit = true;

	Ki67_basic.add_phase_link(0, 1, NULL); // - to +
	Ki67_basic.add_phase_link(1, 0, NULL); // + to -

	Ki67_basic.transition_rate(0, 1) = 1.0 / (4.59 * 60.0); // MCF10A cells are ~4.59 hours in Ki67- state
	Ki67_basic.transition_rate(1, 0) = 1.0 / (15.5 * 60.0); // length of Ki67+ states in advanced model
	Ki67_basic.phase_link(1, 0).fixed_duration = true;

	Ki67_basic.phases[0].entry_function = NULL; // standard_Ki67_negative_phase_entry_function;
	Ki67_basic.phases[1].entry_function = standard_Ki67_positive_phase_entry_function;

	// Ki67_advanced:

	Ki67_advanced.code = constants::advanced_Ki67_cycle_model;
	Ki67_advanced.name = "Ki67 (advanced)";

	Ki67_advanced.data.time_units = "min";

	Ki67_advanced.add_phase(constants::Ki67_negative, "Ki67-");
	Ki67_advanced.add_phase(constants::Ki67_positive_premitotic, "Ki67+ (premitotic)");
	Ki67_advanced.add_phase(constants::Ki67_positive_postmitotic, "Ki67+ (postmitotic)");

	Ki67_advanced.phases[1].division_at_phase_exit = true;

	Ki67_advanced.add_phase_link(0, 1, NULL); // - to +
	Ki67_advanced.add_phase_link(1, 2, NULL); // + (pre-mitotic) to + (post-mitotic)
	Ki67_advanced.add_phase_link(2, 0, NULL); // + to -

	Ki67_advanced.phase_link(1, 2).fixed_duration = true;
	Ki67_advanced.phase_link(2, 0).fixed_duration = true;

	Ki67_advanced.transition_rate(0, 1) = 1.0 / (3.62 * 60.0); // MCF10A cells ~3.62 hours in Ki67- in this fitted model
	Ki67_advanced.transition_rate(1, 2) = 1.0 / (13.0 * 60.0);
	Ki67_advanced.transition_rate(2, 0) = 1.0 / (2.5 * 60.0);

	Ki67_advanced.phases[0].entry_function = NULL; // standard_Ki67_negative_phase_entry_function;
	Ki67_advanced.phases[1].entry_function = standard_Ki67_positive_phase_entry_function;

	return;
}

void create_live_model()
{
	live.code = constants::live_cells_cycle_model;
	live.name = "Live";

	live.data.time_units = "min";

	live.add_phase(constants::live, "Live");

	live.phases[0].division_at_phase_exit = true;

	live.add_phase_link(0, 0, NULL);

	live.transition_rate(0, 0) = 0.0432 / 60.0; // MCF10A have ~0.04 1/hr net birth rate

	live.phases[0].entry_function = standard_live_phase_entry_function;

	return;
}

bool create_cytometry_cycle_models()
{
	// basic one first
	flow_cytometry_cycle_model.code = constants::flow_cytometry_cycle_model;
	flow_cytometry_cycle_model.name = "Flow cytometry model (basic)";

	flow_cytometry_cycle_model.data.time_units = "min";

	flow_cytometry_cycle_model.add_phase(constants::G0G1_phase, "G0/G1");
	flow_cytometry_cycle_model.add_phase(constants::S_phase, "S");
	flow_cytometry_cycle_model.add_phase(constants::G2M_phase, "G2/M");

	flow_cytometry_cycle_model.phases[2].division_at_phase_exit = true;

	flow_cytometry_cycle_model.add_phase_link(0, 1, NULL); // G0/G1 to S
	flow_cytometry_cycle_model.add_phase_link(1, 2, NULL); // S to G2/M
	flow_cytometry_cycle_model.add_phase_link(2, 0, NULL); // G2/M to G0/G1

	// need reference values!
	// https://www.ncbi.nlm.nih.gov/books/NBK9876/
	flow_cytometry_cycle_model.transition_rate(0, 1) = 0.00324; // 5.15 hours in G0/G1 by fitting
	flow_cytometry_cycle_model.transition_rate(1, 2) = 0.00208; // 8 hours in S
	flow_cytometry_cycle_model.transition_rate(2, 0) = 0.00333; // 5 hours in G2/M


	flow_cytometry_cycle_model.phases[0].entry_function = NULL;					  //  ;
	flow_cytometry_cycle_model.phases[1].entry_function = S_phase_entry_function; // Double nuclear volume ;
	flow_cytometry_cycle_model.phases[2].entry_function = NULL;

	// // expanded flow cytometry model

	flow_cytometry_separated_cycle_model.code = constants::flow_cytometry_separated_cycle_model;
	flow_cytometry_separated_cycle_model.name = "Flow cytometry model (separated)";

	flow_cytometry_separated_cycle_model.data.time_units = "min";

	flow_cytometry_separated_cycle_model.add_phase(constants::G0G1_phase, "G0/G1");
	flow_cytometry_separated_cycle_model.add_phase(constants::S_phase, "S");
	flow_cytometry_separated_cycle_model.add_phase(constants::G2_phase, "G2");
	flow_cytometry_separated_cycle_model.add_phase(constants::M_phase, "M");


	flow_cytometry_separated_cycle_model.phases[3].division_at_phase_exit = true;

	flow_cytometry_separated_cycle_model.add_phase_link(0, 1, NULL); // G0/G1 to S
	flow_cytometry_separated_cycle_model.add_phase_link(1, 2, NULL); // S to G2
	flow_cytometry_separated_cycle_model.add_phase_link(2, 3, NULL); // G2 to M
	flow_cytometry_separated_cycle_model.add_phase_link(3, 0, NULL); // M to G0/G1

	// need reference values!
	flow_cytometry_separated_cycle_model.transition_rate(0, 1) = 0.00335; // 4.98 hours in G0/G1
	flow_cytometry_separated_cycle_model.transition_rate(1, 2) = 0.00208; // 8 hours in S
	flow_cytometry_separated_cycle_model.transition_rate(2, 3) = 0.00417; // 4 hours in G2
	flow_cytometry_separated_cycle_model.transition_rate(3, 0) = 0.0167;  // 1 hour in M

	flow_cytometry_separated_cycle_model.phases[0].entry_function = NULL;					//  ;
	flow_cytometry_separated_cycle_model.phases[1].entry_function = S_phase_entry_function; // Double nuclear volume ;
	flow_cytometry_separated_cycle_model.phases[2].entry_function = NULL;
	flow_cytometry_separated_cycle_model.phases[3].entry_function = NULL;

	return true;
}

void create_cycling_quiescent_model()
{
	// Ki67_basic:

	cycling_quiescent.code = constants::cycling_quiescent_model;
	cycling_quiescent.name = "Cycling-Quiescent model";

	cycling_quiescent.data.time_units = "min";

	cycling_quiescent.add_phase(constants::quiescent, "Quiescent");
	cycling_quiescent.add_phase(constants::cycling, "Cycling");

	cycling_quiescent.phases[1].division_at_phase_exit = true;

	cycling_quiescent.add_phase_link(0, 1, NULL); // Q to C
	cycling_quiescent.add_phase_link(1, 0, NULL); // C to Q

	cycling_quiescent.transition_rate(0, 1) = 1.0 / (4.59 * 60.0); // MCF10A cells are ~4.59 hours in Ki67- state
	cycling_quiescent.transition_rate(1, 0) = 1.0 / (15.5 * 60.0); // length of Ki67+ states in advanced model
	cycling_quiescent.phase_link(1, 0).fixed_duration = true;

	cycling_quiescent.phases[0].entry_function = NULL;
	cycling_quiescent.phases[1].entry_function = standard_cycling_entry_function;

	return;
}

bool create_standard_cell_cycle_models()
{
	if (PhysiCell_standard_cycle_models_initialized)
	{
		return false;
	}

	create_ki67_models();
	create_live_model();

	create_cytometry_cycle_models();

	create_cycling_quiescent_model();

	PhysiCell_standard_cycle_models_initialized = true;
	if (PhysiCell_standard_death_models_initialized)
	{
		PhysiCell_standard_models_initialized = true;
	}

	return true;
}

void create_standard_apoptosis_model()
{
	// set default parameters for apoptosis
	apoptosis_parameters.time_units = "min";

	apoptosis_parameters.cytoplasmic_biomass_change_rate = 1.0 / 60.0;
	apoptosis_parameters.nuclear_biomass_change_rate = 0.35 / 60.0;

	apoptosis_parameters.unlysed_fluid_change_rate = 3.0 / 60.0;
	apoptosis_parameters.lysed_fluid_change_rate = 0.0;

	apoptosis_parameters.calcification_rate = 0.0;

	apoptosis_parameters.relative_rupture_volume = 2.0;

	// set up the apoptosis model
	apoptosis.name = "Apoptosis";
	apoptosis.code = constants::apoptosis_death_model;

	// add the main phase for this model, make sure it
	// triggers the appropriate entry function, and note that
	// it should trigger cell removal at its end
	apoptosis.add_phase(constants::apoptotic, "Apoptotic");
	apoptosis.phases[0].entry_function = standard_apoptosis_entry_function;
	apoptosis.phases[0].removal_at_phase_exit = true;

	// add an empty junk debris phase for this model
	apoptosis.add_phase(constants::debris, "Debris");

	// Add a link between these phases. Set the cell to be removed
	// upon this transition. (So the "debris" phase should never be entered).
	apoptosis.add_phase_link(0, 1, NULL);
	apoptosis.transition_rate(0, 1) = 1.0 / (8.6 * 60.0);

	// Use the deterministic model, where this phase has fixed duration
	apoptosis.phase_link(0, 1).fixed_duration = true;

	return;
}

void create_standard_necrosis_model()
{
	// set default parameters for necrosis
	necrosis_parameters.time_units = "min";

	necrosis_parameters.cytoplasmic_biomass_change_rate = 0.0032 / 60.0;
	necrosis_parameters.nuclear_biomass_change_rate = 0.013 / 60.0;

	necrosis_parameters.unlysed_fluid_change_rate = 0.67 / 60.0;
	necrosis_parameters.lysed_fluid_change_rate = 0.050 / 60.0;

	necrosis_parameters.calcification_rate = 0.0042 / 60.0;

	necrosis_parameters.relative_rupture_volume = 2.0;

	// set up the necrosis model
	necrosis.name = "Necrosis";
	necrosis.code = constants::necrosis_death_model;

	necrosis.add_phase(constants::necrotic_swelling, "Necrotic (swelling)");
	necrosis.phases[0].entry_function = standard_necrosis_entry_function;

	necrosis.add_phase(constants::necrotic_lysed, "Necrotic (lysed)");
	necrosis.phases[1].entry_function = standard_lysis_entry_function;
	necrosis.phases[1].removal_at_phase_exit = true;

	// add an empty junk debris phase for this model
	necrosis.add_phase(constants::debris, "Debris");


	necrosis.add_phase_link(0, 1, standard_necrosis_arrest_function);
	necrosis.add_phase_link(1, 2, NULL);

	necrosis.transition_rate(0, 1) = 9e9; // set high so it's always evaluating against the "arrest"
	necrosis.transition_rate(1, 2) = 1.0 / (60.0 * 24.0 * 60.0); // 60 days max

	// Deterministically remove the necrotic cell if it has been 60 days
	necrosis.phase_link(1, 2).fixed_duration = true;

	return;
}

bool create_standard_cell_death_models()
{
	if (PhysiCell_standard_death_models_initialized)
	{
		return false;
	}

	create_standard_apoptosis_model();
	create_standard_necrosis_model();

	PhysiCell_standard_death_models_initialized = true;
	if (PhysiCell_standard_cycle_models_initialized)
	{
		PhysiCell_standard_models_initialized = true;
	}

	return true;
}

bool create_standard_cycle_and_death_models()
{
	bool output = false;
	if (create_standard_cell_cycle_models())
	{
		output = true;
	}
	if (create_standard_cell_death_models())
	{
		output = true;
	}

	return output;
}

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
			cell.convert(*e.cell_definitions[i]);
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
		standard_cell_transformations(*cell, e);

		// New March 2023 in Version 1.12.0
		// call the rules-based code to update the phenotype
		if (e.rules_enabled)
		{
			apply_ruleset(cell.get(), e);
		}
		if (get_single_signal(cell.get(), "necrotic", e) > 0.5)
		{
			real_t rupture = cell->phenotype.volume.rupture_volume();
			real_t volume = cell->phenotype.volume.total();
			if (volume > rupture)
			{
				std::cout << cell->phenotype.volume.total() << " vs " << cell->phenotype.volume.rupture_volume()
						  << " dead: " << get_single_signal(cell.get(), "dead", e) << std::endl;
				std::cout << cell->phenotype.cycle.current_phase_index() << " "
						  << cell->phenotype.cycle.pCycle_Model->name << std::endl;
			}
		}

		// call the custom code to update the phenotype
		if (cell->functions.update_phenotype)
		{
			cell->functions.update_phenotype(*cell);
		}

		// update volume
		if (cell->functions.volume_update_function)
		{
			cell->functions.volume_update_function(*cell);
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

template void advanced_chemotaxis_function<true>(cell& cell);
template void advanced_chemotaxis_function<false>(cell& cell);

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

void standard_elastic_contact_function(cell& lhs, cell& rhs)
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

void update_cell_and_death_parameters_O2_based(cell& cell)
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
	static index_t start_phase_index; // Q_phase_index;
	static index_t end_phase_index;	  // K_phase_index;
	static index_t necrosis_index;

	static index_t oxygen_substrate_index = cell.e().m.find_substrate_index("oxygen");

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

	real_t pO2 = (cell.nearest_density_vector())[oxygen_substrate_index]; // constants::oxygen_index];

	// this multiplier is for linear interpolation of the oxygen value
	real_t multiplier = 1.0;
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

void initialize_default_cell_definition(environment& e)
{
	// If the standard models have not yet been created, do so now.
	create_standard_cycle_and_death_models();

	e.create_cell_definition();

	e.cell_defaults().type = 0;
	e.cell_defaults().name = "breast epithelium";

	// set up the default functions
	e.cell_defaults().phenotype.cycle.sync_to_cycle_model(Ki67_advanced);

	e.cell_defaults().functions.volume_update_function = standard_volume_update_function;

	e.cell_defaults().functions.update_phenotype = update_cell_and_death_parameters_O2_based; // NULL;

	// add the standard death models to the default phenotype.
	e.cell_defaults().phenotype.death.add_death_model(0.00319 / 60.0, &apoptosis, apoptosis_parameters);
	// MCF10A, to get a 2% apoptotic index
	e.cell_defaults().phenotype.death.add_death_model(0.0, &necrosis, necrosis_parameters);

	// Cell_Parameters, Custom_Cell_Data, Cell_Functions are handled by the default constructor

	e.cell_defaults().phenotype.volume.set_defaults();
	e.cell_defaults().phenotype.geometry.set_defaults();
	e.cell_defaults().phenotype.mechanics.set_defaults();
	e.cell_defaults().phenotype.motility.set_defaults();
	e.cell_defaults().phenotype.secretion.set_defaults();
	e.cell_defaults().phenotype.molecular.set_defaults();
	e.cell_defaults().phenotype.cell_interactions.set_defaults();
	e.cell_defaults().phenotype.cell_transformations.set_defaults();
}

} // namespace physicell
