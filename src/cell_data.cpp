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

void volume_data::add(index_t size)
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
}

void geometry_data::add(index_t size)
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
}

void mechanics_data::add(index_t size, index_t cell_definitions_count)
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
}

void motility_data::add(index_t index, index_t dims, index_t substrates_count)
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
}

cell_data::cell_data(environment& e) : agent_data(e.m), agents_count(agent_data.agents_count), e(e) {}

void cell_data::add()
{
	agent_data.add();

	volumes.add(agents_count);
	geometries.add(agents_count);
	mechanics.add(agents_count, e.cell_definitions_count);
	motility.add(agents_count, e.mechanics_mesh.dims, e.m.substrates_count);

	neighbors.resize(agents_count);
	springs.resize(agents_count);

	previous_velocities.resize(agents_count * e.m.mesh.dims, 0);
	velocities.resize(agents_count * e.m.mesh.dims, 0);
	cell_definition_indices.resize(agents_count, 0);
	simple_pressures.resize(agents_count, 0);
	is_movable.resize(agents_count, 0);
}

void cell_data::remove(index_t index)
{
	agent_data.remove(index);

	if (index == agents_count)
	{
		return;
	}

	volumes.remove(index, agents_count);
	geometries.remove(index, agents_count);
	mechanics.remove(index, agents_count, e.cell_definitions_count);
	motility.remove(index, agents_count, e.mechanics_mesh.dims, e.m.substrates_count);

	neighbors[index] = std::move(neighbors[agents_count]);
	springs[index] = std::move(springs[agents_count]);

	move_vector(previous_velocities.data() + index * e.m.mesh.dims,
				previous_velocities.data() + agents_count * e.m.mesh.dims, e.m.mesh.dims);
	move_vector(velocities.data() + index * e.m.mesh.dims, velocities.data() + agents_count * e.m.mesh.dims,
				e.m.mesh.dims);
	move_scalar(cell_definition_indices.data() + index, cell_definition_indices.data() + agents_count);
	move_scalar(simple_pressures.data() + index, simple_pressures.data() + agents_count);
	move_scalar(is_movable.data() + index, is_movable.data() + agents_count);
}
