#include "cell_data.h"

#include <cstring>

#include <BioFVM/microenvironment.h>

#include "environment.h"

using namespace biofvm;
using namespace physicell;

template <typename T>
void move_scalar(T* dst, const T* src)
{
	dst[0] = src[0];
}

template <typename T>
void move_vector(T* dst, const T* src, biofvm::index_t size)
{
	std::memcpy(dst, src, size * sizeof(T));
}

void volume_data::resize(index_t size)
{
	solid.resize(size, 0);
	fluid.resize(size, 0);
	fluid_fraction.resize(size, 0);

	nuclear.resize(size, 0);
	nuclear_fluid.resize(size, 0);
	nuclear_solid.resize(size, 0);

	cytoplasmic.resize(size, 0);
	cytoplasmic_fluid.resize(size, 0);
	cytoplasmic_solid.resize(size, 0);

	calcified_fraction.resize(size, 0);

	cytoplasmic_to_nuclear_ratio.resize(size, 0);

	rupture_volume.resize(size, 0);

	cytoplasmic_biomass_change_rate.resize(size, 0);
	nuclear_biomass_change_rate.resize(size, 0);
	fluid_change_rate.resize(size, 0);
	calcification_rate.resize(size, 0);
	target_solid_cytoplasmic.resize(size, 0);
	target_solid_nuclear.resize(size, 0);
	target_fluid_fraction.resize(size, 0);
	target_cytoplasmic_to_nuclear_ratio.resize(size, 0);
	relative_rupture_volume.resize(size, 0);
}

void volume_data::remove(index_t index, index_t size)
{
	move_scalar(solid.data() + index, solid.data() + size);
	move_scalar(fluid.data() + index, fluid.data() + size);
	move_scalar(fluid_fraction.data() + index, fluid_fraction.data() + size);

	move_scalar(nuclear.data() + index, nuclear.data() + size);
	move_scalar(nuclear_fluid.data() + index, nuclear_fluid.data() + size);
	move_scalar(nuclear_solid.data() + index, nuclear_solid.data() + size);

	move_scalar(cytoplasmic.data() + index, cytoplasmic.data() + size);
	move_scalar(cytoplasmic_fluid.data() + index, cytoplasmic_fluid.data() + size);
	move_scalar(cytoplasmic_solid.data() + index, cytoplasmic_solid.data() + size);

	move_scalar(calcified_fraction.data() + index, calcified_fraction.data() + size);

	move_scalar(cytoplasmic_to_nuclear_ratio.data() + index, cytoplasmic_to_nuclear_ratio.data() + size);

	move_scalar(rupture_volume.data() + index, rupture_volume.data() + size);

	move_scalar(cytoplasmic_biomass_change_rate.data() + index, cytoplasmic_biomass_change_rate.data() + size);
	move_scalar(nuclear_biomass_change_rate.data() + index, nuclear_biomass_change_rate.data() + size);
	move_scalar(fluid_change_rate.data() + index, fluid_change_rate.data() + size);
	move_scalar(calcification_rate.data() + index, calcification_rate.data() + size);
	move_scalar(target_solid_cytoplasmic.data() + index, target_solid_cytoplasmic.data() + size);
	move_scalar(target_solid_nuclear.data() + index, target_solid_nuclear.data() + size);
	move_scalar(target_fluid_fraction.data() + index, target_fluid_fraction.data() + size);
	move_scalar(target_cytoplasmic_to_nuclear_ratio.data() + index, target_cytoplasmic_to_nuclear_ratio.data() + size);
	move_scalar(relative_rupture_volume.data() + index, relative_rupture_volume.data() + size);

	resize(size - 1);
}

void geometry_data::resize(index_t size)
{
	radius.resize(size, 0);
	nuclear_radius.resize(size, 0);
	surface_area.resize(size, 0);
	polarity.resize(size, 0);
}

void geometry_data::remove(index_t index, index_t size)
{
	move_scalar(radius.data() + index, radius.data() + size);
	move_scalar(nuclear_radius.data() + index, nuclear_radius.data() + size);
	move_scalar(surface_area.data() + index, surface_area.data() + size);
	move_scalar(polarity.data() + index, polarity.data() + size);

	resize(size - 1);
}

void mechanics_data::resize(index_t size, index_t cell_definitions_count)
{
	cell_cell_adhesion_strength.resize(size, 0);
	cell_BM_adhesion_strength.resize(size, 0);

	cell_cell_repulsion_strength.resize(size, 0);
	cell_BM_repulsion_strength.resize(size, 0);

	cell_adhesion_affinities.resize(size * cell_definitions_count, 0);

	relative_maximum_adhesion_distance.resize(size, 0);

	maximum_number_of_attachments.resize(size, 0);
	attachment_elastic_constant.resize(size, 0);

	attachment_rate.resize(size, 0);
	detachment_rate.resize(size, 0);
}

void mechanics_data::remove(index_t index, index_t size, index_t cell_definitions_count)
{
	move_scalar(cell_cell_adhesion_strength.data() + index, cell_cell_adhesion_strength.data() + size);
	move_scalar(cell_BM_adhesion_strength.data() + index, cell_BM_adhesion_strength.data() + size);

	move_scalar(cell_cell_repulsion_strength.data() + index, cell_cell_repulsion_strength.data() + size);
	move_scalar(cell_BM_repulsion_strength.data() + index, cell_BM_repulsion_strength.data() + size);

	move_vector(cell_adhesion_affinities.data() + index * cell_definitions_count,
				cell_adhesion_affinities.data() + size * cell_definitions_count, cell_definitions_count);

	move_scalar(relative_maximum_adhesion_distance.data() + index, relative_maximum_adhesion_distance.data() + size);

	move_scalar(maximum_number_of_attachments.data() + index, maximum_number_of_attachments.data() + size);
	move_scalar(attachment_elastic_constant.data() + index, attachment_elastic_constant.data() + size);

	move_scalar(attachment_rate.data() + index, attachment_rate.data() + size);
	move_scalar(detachment_rate.data() + index, detachment_rate.data() + size);

	resize(size - 1, cell_definitions_count);
}

void motility_data::resize(index_t index, index_t dims, index_t substrates_count)
{
	is_motile.resize(index, 0);
	persistence_time.resize(index, 0);
	migration_speed.resize(index, 0);

	migration_bias_direction.resize(index * dims, 0);
	migration_bias.resize(index, 0);

	motility_vector.resize(index * dims, 0);

	restrict_to_2d.resize(index, 0);

	chemotaxis_index.resize(index, 0);
	chemotaxis_direction.resize(index, 0);
	chemotactic_sensitivities.resize(index * substrates_count, 0);

	update_migration_bias_direction.resize(index, nullptr);
}

void motility_data::remove(index_t index, index_t size, index_t dims, index_t substrates_count)
{
	move_scalar(is_motile.data() + index, is_motile.data() + size);
	move_scalar(persistence_time.data() + index, persistence_time.data() + size);
	move_scalar(migration_speed.data() + index, migration_speed.data() + size);

	move_vector(migration_bias_direction.data() + index * dims, migration_bias_direction.data() + size * dims, dims);
	move_scalar(migration_bias.data() + index, migration_bias.data() + size);

	move_vector(motility_vector.data() + index * dims, motility_vector.data() + size * dims, dims);

	move_scalar(restrict_to_2d.data() + index, restrict_to_2d.data() + size);

	move_scalar(chemotaxis_index.data() + index, chemotaxis_index.data() + size);
	move_scalar(chemotaxis_direction.data() + index, chemotaxis_direction.data() + size);
	move_vector(chemotactic_sensitivities.data() + index * substrates_count,
				chemotactic_sensitivities.data() + size * substrates_count, substrates_count);

	move_scalar(update_migration_bias_direction.data() + index, update_migration_bias_direction.data() + size);

	resize(size - 1, dims, substrates_count);
}

void interactions_data::resize(index_t size, index_t cell_definitions_count)
{
	dead_phagocytosis_rate.resize(size, 0);
	live_phagocytosis_rates.resize(size * cell_definitions_count, 0);

	damage_rate.resize(size, 0);
	attack_rates.resize(size * cell_definitions_count, 0);
	immunogenicities.resize(size * cell_definitions_count, 0);

	fusion_rates.resize(size * cell_definitions_count, 0);
}

void interactions_data::remove(index_t index, index_t size, index_t cell_definitions_count)
{
	move_scalar(dead_phagocytosis_rate.data() + index, dead_phagocytosis_rate.data() + size);
	move_vector(live_phagocytosis_rates.data() + index * cell_definitions_count,
				live_phagocytosis_rates.data() + size * cell_definitions_count, cell_definitions_count);

	move_scalar(damage_rate.data() + index, damage_rate.data() + size);
	move_vector(attack_rates.data() + index * cell_definitions_count,
				attack_rates.data() + size * cell_definitions_count, cell_definitions_count);
	move_vector(immunogenicities.data() + index * cell_definitions_count,
				immunogenicities.data() + size * cell_definitions_count, cell_definitions_count);

	move_vector(fusion_rates.data() + index * cell_definitions_count,
				fusion_rates.data() + size * cell_definitions_count, cell_definitions_count);

	resize(size - 1, cell_definitions_count);
}

void death_data::resize(index_t size) { dead.resize(size, 0); }

void death_data::remove(index_t index, index_t size)
{
	move_scalar(dead.data() + index, dead.data() + size);

	resize(size - 1);
}

void transformations_data::resize(index_t size, index_t cell_definitions_count)
{
	transformation_rates.resize(size * cell_definitions_count, 0);
}

void transformations_data::remove(index_t index, index_t size, index_t cell_definitions_count)
{
	move_vector(transformation_rates.data() + index * cell_definitions_count,
				transformation_rates.data() + size * cell_definitions_count, cell_definitions_count);

	resize(size - 1, cell_definitions_count);
}

void cell_state_data::resize(index_t size, index_t dims)
{
	neighbors.resize(size);

	springs.resize(size);
	attached_cells.resize(size);

	orientation.resize(size * dims, 0);
	simple_pressure.resize(size, 0);
	number_of_nuclei.resize(size, 0);

	damage.resize(size, 0);
	total_attack_time.resize(size, 0);
}

void cell_state_data::remove(index_t index, index_t size, index_t dims)
{
	neighbors[index] = std::move(neighbors[size]);
	springs[index] = std::move(springs[size]);
	attached_cells[index] = std::move(attached_cells[size]);

	move_vector(orientation.data() + index * dims, orientation.data() + size * dims, dims);
	move_scalar(simple_pressure.data() + index, simple_pressure.data() + size);
	move_scalar(number_of_nuclei.data() + index, number_of_nuclei.data() + size);

	move_scalar(damage.data() + index, damage.data() + size);
	move_scalar(total_attack_time.data() + index, total_attack_time.data() + size);

	resize(size - 1, dims);
}

cell_data::cell_data(environment& e) : agent_data(e.m), agents_count(agent_data.agents_count), e(e) {}

void cell_data::add()
{
	agent_data.add();
	
	resize();
}

void cell_data::resize()
{
	volumes.resize(agents_count);
	geometries.resize(agents_count);
	mechanics.resize(agents_count, e.cell_definitions_count);
	motilities.resize(agents_count, e.mechanics_mesh.dims, e.m.substrates_count);
	deaths.resize(agents_count);

	states.resize(agents_count, e.m.mesh.dims);

	interactions.resize(agents_count, e.cell_definitions_count);
	transformations.resize(agents_count, e.cell_definitions_count);

	previous_velocities.resize(agents_count * e.m.mesh.dims, 0);
	velocities.resize(agents_count * e.m.mesh.dims, 0);
	cell_definition_indices.resize(agents_count, 0);
	is_movable.resize(agents_count, 0);

	flags.resize(agents_count, cell_state_flag::none);
}

void cell_data::remove(index_t index)
{
	agent_data.remove(index);

	volumes.remove(index, agents_count);
	geometries.remove(index, agents_count);
	mechanics.remove(index, agents_count, e.cell_definitions_count);
	motilities.remove(index, agents_count, e.mechanics_mesh.dims, e.m.substrates_count);
	deaths.remove(index, agents_count);

	states.remove(index, agents_count, e.mechanics_mesh.dims);

	interactions.remove(index, agents_count, e.cell_definitions_count);
	transformations.remove(index, agents_count, e.cell_definitions_count);

	move_vector(previous_velocities.data() + index * e.m.mesh.dims,
				previous_velocities.data() + agents_count * e.m.mesh.dims, e.m.mesh.dims);
	move_vector(velocities.data() + index * e.m.mesh.dims, velocities.data() + agents_count * e.m.mesh.dims,
				e.m.mesh.dims);
	move_scalar(cell_definition_indices.data() + index, cell_definition_indices.data() + agents_count);
	move_scalar(is_movable.data() + index, is_movable.data() + agents_count);

	move_scalar(flags.data() + index, flags.data() + agents_count);
}
