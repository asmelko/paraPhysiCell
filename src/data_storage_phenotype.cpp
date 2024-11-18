#include "data_storage_phenotype.h"

#include "cell.h"
#include "cell_data.h"
#include "environment.h"
#include "solver/host/interactions_solver.h"

using namespace biofvm;
using namespace physicell;

phenotype_data_storage::phenotype_data_storage(cell_data& data, const index_t& index) : data_(data), index_(index) {}

volume_t::volume_t(cell_data& data, const index_t& index) : phenotype_data_storage(data, index) {}

real_t& volume_t::total() { return data_.agent_data.volumes[index_]; }

real_t& volume_t::solid() { return data_.volumes.solid[index_]; }

real_t& volume_t::fluid() { return data_.volumes.fluid[index_]; }

real_t& volume_t::fluid_fraction() { return data_.volumes.fluid_fraction[index_]; }

real_t& volume_t::nuclear() { return data_.volumes.nuclear[index_]; }

real_t& volume_t::nuclear_fluid() { return data_.volumes.nuclear_fluid[index_]; }

real_t& volume_t::nuclear_solid() { return data_.volumes.nuclear_solid[index_]; }

real_t& volume_t::cytoplasmic() { return data_.volumes.cytoplasmic[index_]; }

real_t& volume_t::cytoplasmic_fluid() { return data_.volumes.cytoplasmic_fluid[index_]; }

real_t& volume_t::cytoplasmic_solid() { return data_.volumes.cytoplasmic_solid[index_]; }

real_t& volume_t::calcified_fraction() { return data_.volumes.calcified_fraction[index_]; }

real_t& volume_t::cytoplasmic_to_nuclear_ratio() { return data_.volumes.cytoplasmic_to_nuclear_ratio[index_]; }

real_t& volume_t::rupture_volume() { return data_.volumes.rupture_volume[index_]; }

real_t& volume_t::cytoplasmic_biomass_change_rate() { return data_.volumes.cytoplasmic_biomass_change_rate[index_]; }

real_t& volume_t::nuclear_biomass_change_rate() { return data_.volumes.nuclear_biomass_change_rate[index_]; }

real_t& volume_t::fluid_change_rate() { return data_.volumes.fluid_change_rate[index_]; }

real_t& volume_t::calcification_rate() { return data_.volumes.calcification_rate[index_]; }

real_t& volume_t::target_solid_cytoplasmic() { return data_.volumes.target_solid_cytoplasmic[index_]; }

real_t& volume_t::target_solid_nuclear() { return data_.volumes.target_solid_nuclear[index_]; }

real_t& volume_t::target_fluid_fraction() { return data_.volumes.target_fluid_fraction[index_]; }

real_t& volume_t::target_cytoplasmic_to_nuclear_ratio()
{
	return data_.volumes.target_cytoplasmic_to_nuclear_ratio[index_];
}

real_t& volume_t::relative_rupture_volume() { return data_.volumes.relative_rupture_volume[index_]; }

void volume_t::set_defaults()
{
	// reference parameter values for MCF-7, in cubic microns
	fluid_fraction() = 0.75;

	total() = 2494;
	fluid() = fluid_fraction() * total();
	solid() = total() - fluid();

	nuclear() = 540.0;


	nuclear_fluid() = fluid_fraction() * nuclear();
	nuclear_solid() = nuclear() - nuclear_fluid();

	cytoplasmic() = total() - nuclear();
	cytoplasmic_fluid() = fluid_fraction() * cytoplasmic();
	cytoplasmic_solid() = cytoplasmic() - cytoplasmic_fluid();

	// rates are in units of 1/min
	cytoplasmic_biomass_change_rate() = 0.27 / 60.0;
	nuclear_biomass_change_rate() = 0.33 / 60.0;
	fluid_change_rate() = 3.0 / 60.0;

	calcified_fraction() = 0.0;

	calcification_rate() = 0.0;

	target_solid_cytoplasmic() = cytoplasmic_solid();
	target_solid_nuclear() = nuclear_solid();
	target_fluid_fraction() = fluid_fraction();

	cytoplasmic_to_nuclear_ratio() = cytoplasmic() / (1e-16 + nuclear());
	target_cytoplasmic_to_nuclear_ratio() = cytoplasmic_to_nuclear_ratio();

	// the cell bursts at these volumes
	relative_rupture_volume() = 2.0;
	// as fraction of volume at entry to the current phase
	rupture_volume() = relative_rupture_volume() * total(); // in volume units
}

geometry_t::geometry_t(cell_data& data, const index_t& index) : phenotype_data_storage(data, index) {}

real_t& geometry_t::radius() { return data_.geometries.radius[index_]; }

real_t& geometry_t::nuclear_radius() { return data_.geometries.nuclear_radius[index_]; }

real_t& geometry_t::surface_area() { return data_.geometries.surface_area[index_]; }

real_t& geometry_t::polarity() { return data_.geometries.polarity[index_]; }

void geometry_t::update()
{
	update_geometry(index_, data_.geometries.radius.data(), data_.geometries.nuclear_radius.data(),
					data_.geometries.surface_area.data(), data_.agent_data.volumes.data(),
					data_.volumes.nuclear.data());
}

void geometry_t::set_defaults()
{
	radius() = 8.412710547954228;
	nuclear_radius() = 5.051670902881889;
	surface_area() = 889.3685284131693;

	polarity() = 0.0;
}

mechanics_t::mechanics_t(cell_data& data, const index_t& index) : phenotype_data_storage(data, index) {}

real_t& mechanics_t::cell_cell_adhesion_strength() { return data_.mechanics.cell_cell_adhesion_strength[index_]; }

real_t& mechanics_t::cell_BM_adhesion_strength() { return data_.mechanics.cell_BM_adhesion_strength[index_]; }

real_t& mechanics_t::cell_cell_repulsion_strength() { return data_.mechanics.cell_cell_repulsion_strength[index_]; }

real_t& mechanics_t::cell_BM_repulsion_strength() { return data_.mechanics.cell_BM_repulsion_strength[index_]; }

real_t* mechanics_t::cell_adhesion_affinities()
{
	return data_.mechanics.cell_adhesion_affinities.data() + index_ * data_.e.cell_definitions_count;
}

real_t& mechanics_t::relative_maximum_adhesion_distance()
{
	return data_.mechanics.relative_maximum_adhesion_distance[index_];
}

index_t& mechanics_t::maximum_number_of_attachments() { return data_.mechanics.maximum_number_of_attachments[index_]; }

real_t& mechanics_t::attachment_elastic_constant() { return data_.mechanics.attachment_elastic_constant[index_]; }

real_t& mechanics_t::attachment_rate() { return data_.mechanics.attachment_rate[index_]; }

real_t& mechanics_t::detachment_rate() { return data_.mechanics.detachment_rate[index_]; }

void mechanics_t::set_relative_equilibrium_distance(real_t new_value)
{
	if (new_value > 2.0)
	{
		std::cout << "**** Warning in function " << __FUNCTION__ << " in " << __FILE__ << " : " << std::endl
				  << "\tAttempted to set equilibrium distance exceeding two cell radii." << std::endl
				  << "\tWe will cap the equilibrium distance at 2.0 cell radii." << std::endl
				  << "****" << std::endl
				  << std::endl;

		new_value = 2.0;
	}

	// adjust the adhesive coefficient to achieve the new (relative) equilibrium distance

	real_t temp1 = new_value;
	temp1 /= 2.0;

	real_t temp2 = 1.0;
	temp2 -= temp1; // 1 - s_relative/2.0

	temp1 /= relative_maximum_adhesion_distance(); // s_relative/(2*relative_maximum_adhesion_distance);
	temp1 *= -1.0;								   // -s_relative/(2*relative_maximum_adhesion_distance);
	temp1 += 1.0;								   // 1.0 -s_relative/(2*relative_maximum_adhesion_distance);

	temp2 /= temp1;
	temp2 *= temp2;

	cell_cell_adhesion_strength() = cell_cell_repulsion_strength();
	cell_cell_adhesion_strength() *= temp2;
}

void mechanics_t::set_absolute_equilibrium_distance(real_t radius, real_t new_value)
{
	set_relative_equilibrium_distance(new_value / radius);
}

void mechanics_t::set_defaults()
{
	cell_cell_adhesion_strength() = 0.4;
	cell_BM_adhesion_strength() = 4.0;

	cell_cell_repulsion_strength() = 10.0;
	cell_BM_repulsion_strength() = 100.0;

	std::fill(cell_adhesion_affinities(), cell_adhesion_affinities() + data_.e.cell_definitions_count, 1);

	// this is a multiple of the cell (equivalent) radius
	relative_maximum_adhesion_distance() = 1.25;
	// maximum_adhesion_distance = 0.0;

	/* for spring attachments */
	maximum_number_of_attachments() = 12;
	attachment_elastic_constant() = 0.01;

	attachment_rate() = 0; // 10.0 prior ot March 2023
	detachment_rate() = 0;
}

motility_t::motility_t(cell_data& data, const index_t& index) : phenotype_data_storage(data, index) {}

std::uint8_t& motility_t::is_motile() { return data_.motilities.is_motile[index_]; }

real_t& motility_t::persistence_time() { return data_.motilities.persistence_time[index_]; }

real_t& motility_t::migration_speed() { return data_.motilities.migration_speed[index_]; }

real_t* motility_t::migration_bias_direction()
{
	return data_.motilities.migration_bias_direction.data() + index_ * data_.e.mechanics_mesh.dims;
}

real_t& motility_t::migration_bias() { return data_.motilities.migration_bias[index_]; }

real_t* motility_t::motility_vector()
{
	return data_.motilities.motility_vector.data() + index_ * data_.e.mechanics_mesh.dims;
}

index_t& motility_t::chemotaxis_index() { return data_.motilities.chemotaxis_index[index_]; }

index_t& motility_t::chemotaxis_direction() { return data_.motilities.chemotaxis_direction[index_]; }

real_t* motility_t::chemotactic_sensitivities()
{
	return data_.motilities.chemotactic_sensitivities.data() + index_ * data_.e.m.substrates_count;
}

void motility_t::set_defaults()
{
	is_motile() = false;

	persistence_time() = 1.0;
	migration_speed() = 1.0;

	std::fill(migration_bias_direction(), migration_bias_direction() + data_.e.mechanics_mesh.dims, 0);
	migration_bias() = 0.0;

	update_migration_bias_direction() = nullptr;

	std::fill(motility_vector(), motility_vector() + data_.e.mechanics_mesh.dims, 0);

	chemotaxis_index() = 0;
	chemotaxis_direction() = 1;
	std::fill(chemotactic_sensitivities(), chemotactic_sensitivities() + data_.e.m.substrates_count, 0);
}

motility_data::direction_update_func& motility_t::update_migration_bias_direction()
{
	return data_.motilities.update_migration_bias_direction[index_];
}

molecular_t::molecular_t(cell_data& data, const index_t& index) : phenotype_data_storage(data, index) {}

real_t* molecular_t::internalized_total_substrates()
{
	return data_.agent_data.internalized_substrates.data() + index_ * data_.e.m.substrates_count;
}

real_t* molecular_t::fraction_released_at_death()
{
	return data_.agent_data.fraction_released_at_death.data() + index_ * data_.e.m.substrates_count;
}

real_t* molecular_t::fraction_transferred_when_ingested()
{
	return data_.agent_data.fraction_transferred_when_ingested.data() + index_ * data_.e.m.substrates_count;
}

void molecular_t::set_defaults()
{
	std::fill(internalized_total_substrates(), internalized_total_substrates() + data_.e.m.substrates_count, 0);
	std::fill(fraction_released_at_death(), fraction_released_at_death() + data_.e.m.substrates_count, 0);
	std::fill(fraction_transferred_when_ingested(), fraction_transferred_when_ingested() + data_.e.m.substrates_count,
			  0);
}

interactions_t::interactions_t(cell_data& data, const index_t& index) : phenotype_data_storage(data, index) {}

real_t& interactions_t::apoptotic_phagocytosis_rate() { return data_.interactions.apoptotic_phagocytosis_rate[index_]; }

real_t& interactions_t::necrotic_phagocytosis_rate() { return data_.interactions.necrotic_phagocytosis_rate[index_]; }

real_t& interactions_t::other_dead_phagocytosis_rate()
{
	return data_.interactions.other_dead_phagocytosis_rate[index_];
}

real_t* interactions_t::live_phagocytosis_rates()
{
	return data_.interactions.live_phagocytosis_rates.data() + index_ * data_.e.cell_definitions_count;
}

real_t* interactions_t::attack_rates()
{
	return data_.interactions.attack_rates.data() + index_ * data_.e.cell_definitions_count;
}

real_t& interactions_t::attack_damage_rate() { return data_.interactions.attack_damage_rates[index_]; }

real_t& interactions_t::attack_duration() { return data_.interactions.attack_durations[index_]; }

real_t& interactions_t::total_damage_delivered() { return data_.interactions.total_damage_delivered[index_]; }

index_t& interactions_t::attack_target() { return data_.interactions.attack_targets[index_]; }

real_t* interactions_t::immunogenicities()
{
	return data_.interactions.immunogenicities.data() + index_ * data_.e.cell_definitions_count;
}

real_t* interactions_t::fusion_rates()
{
	return data_.interactions.fusion_rates.data() + index_ * data_.e.cell_definitions_count;
}

void interactions_t::set_defaults()
{
	apoptotic_phagocytosis_rate() = 0.0;
	necrotic_phagocytosis_rate() = 0.0;
	other_dead_phagocytosis_rate() = 0.0;
	std::fill(live_phagocytosis_rates(), live_phagocytosis_rates() + data_.e.cell_definitions_count, 0);
	std::fill(attack_rates(), attack_rates() + data_.e.cell_definitions_count, 0);
	attack_damage_rate() = 1.0;
	attack_duration() = 30.0;
	total_damage_delivered() = 0.0;
	attack_target() = no_target;
	std::fill(immunogenicities(), immunogenicities() + data_.e.cell_definitions_count, 1);
	std::fill(fusion_rates(), fusion_rates() + data_.e.cell_definitions_count, 0);
}

secretion_t::secretion_t(cell_data& data, const index_t& index) : phenotype_data_storage(data, index) {}

real_t* secretion_t::secretion_rates()
{
	return data_.agent_data.secretion_rates.data() + index_ * data_.e.m.substrates_count;
}

real_t* secretion_t::uptake_rates()
{
	return data_.agent_data.uptake_rates.data() + index_ * data_.e.m.substrates_count;
}

real_t* secretion_t::saturation_densities()
{
	return data_.agent_data.saturation_densities.data() + index_ * data_.e.m.substrates_count;
}

real_t* secretion_t::net_export_rates()
{
	return data_.agent_data.net_export_rates.data() + index_ * data_.e.m.substrates_count;
}

void secretion_t::set_all_secretion_to_zero()
{
	for (index_t i = 0; i < data_.e.m.substrates_count; i++)
	{
		secretion_rates()[i] = 0;
	}
}

void secretion_t::scale_all_uptake_by_factor(real_t factor)
{
	for (index_t i = 0; i < data_.e.m.substrates_count; i++)
	{
		uptake_rates()[i] *= factor;
	}
}

void secretion_t::set_defaults()
{
	std::fill(secretion_rates(), secretion_rates() + data_.e.m.substrates_count, 0);
	std::fill(uptake_rates(), uptake_rates() + data_.e.m.substrates_count, 0);
	std::fill(saturation_densities(), saturation_densities() + data_.e.m.substrates_count, 0);
	std::fill(net_export_rates(), net_export_rates() + data_.e.m.substrates_count, 0);
}

transformations_t::transformations_t(cell_data& data, const index_t& index) : phenotype_data_storage(data, index) {}

real_t* transformations_t::transformation_rates()
{
	return data_.transformations.transformation_rates.data() + index_ * data_.e.cell_definitions_count;
}

void transformations_t::set_defaults()
{
	std::fill(transformation_rates(), transformation_rates() + data_.e.cell_definitions_count, 0);
}

void volume_t::copy(volume_t& dest)
{
	dest.total() = total();
	dest.solid() = solid();
	dest.fluid() = fluid();
	dest.fluid_fraction() = fluid_fraction();
	dest.nuclear() = nuclear();
	dest.nuclear_fluid() = nuclear_fluid();
	dest.nuclear_solid() = nuclear_solid();
	dest.cytoplasmic() = cytoplasmic();
	dest.cytoplasmic_fluid() = cytoplasmic_fluid();
	dest.cytoplasmic_solid() = cytoplasmic_solid();
	dest.calcified_fraction() = calcified_fraction();
	dest.cytoplasmic_to_nuclear_ratio() = cytoplasmic_to_nuclear_ratio();
	dest.rupture_volume() = rupture_volume();

	dest.cytoplasmic_biomass_change_rate() = cytoplasmic_biomass_change_rate();
	dest.nuclear_biomass_change_rate() = nuclear_biomass_change_rate();
	dest.fluid_change_rate() = fluid_change_rate();
	dest.calcification_rate() = calcification_rate();
	dest.target_solid_cytoplasmic() = target_solid_cytoplasmic();
	dest.target_solid_nuclear() = target_solid_nuclear();
	dest.target_fluid_fraction() = target_fluid_fraction();
	dest.target_cytoplasmic_to_nuclear_ratio() = target_cytoplasmic_to_nuclear_ratio();
	dest.relative_rupture_volume() = relative_rupture_volume();
}

void geometry_t::copy(geometry_t& dest)
{
	dest.radius() = radius();
	dest.nuclear_radius() = nuclear_radius();
	dest.surface_area() = surface_area();
	dest.polarity() = polarity();
}

void mechanics_t::copy(mechanics_t& dest)
{
	dest.cell_cell_adhesion_strength() = cell_cell_adhesion_strength();
	dest.cell_BM_adhesion_strength() = cell_BM_adhesion_strength();
	dest.cell_cell_repulsion_strength() = cell_cell_repulsion_strength();
	dest.cell_BM_repulsion_strength() = cell_BM_repulsion_strength();
	for (index_t i = 0; i < data_.e.cell_definitions_count; i++)
	{
		dest.cell_adhesion_affinities()[i] = cell_adhesion_affinities()[i];
	}
	dest.relative_maximum_adhesion_distance() = relative_maximum_adhesion_distance();
	dest.maximum_number_of_attachments() = maximum_number_of_attachments();
	dest.attachment_elastic_constant() = attachment_elastic_constant();
	dest.attachment_rate() = attachment_rate();
	dest.detachment_rate() = detachment_rate();
}

void motility_t::copy(motility_t& dest)
{
	dest.is_motile() = is_motile();
	dest.persistence_time() = persistence_time();
	dest.migration_speed() = migration_speed();
	for (index_t i = 0; i < data_.e.mechanics_mesh.dims; i++)
	{
		dest.migration_bias_direction()[i] = migration_bias_direction()[i];
	}
	dest.migration_bias() = migration_bias();
	for (index_t i = 0; i < data_.e.mechanics_mesh.dims; i++)
	{
		dest.motility_vector()[i] = motility_vector()[i];
	}
	dest.chemotaxis_index() = chemotaxis_index();
	dest.chemotaxis_direction() = chemotaxis_direction();
	for (index_t i = 0; i < data_.e.m.substrates_count; i++)
	{
		dest.chemotactic_sensitivities()[i] = chemotactic_sensitivities()[i];
	}

	dest.update_migration_bias_direction() = update_migration_bias_direction();
}

void molecular_t::copy(molecular_t& dest)
{
	for (index_t i = 0; i < data_.e.m.substrates_count; i++)
	{
		dest.internalized_total_substrates()[i] = internalized_total_substrates()[i];
	}
	for (index_t i = 0; i < data_.e.m.substrates_count; i++)
	{
		dest.fraction_released_at_death()[i] = fraction_released_at_death()[i];
	}
	for (index_t i = 0; i < data_.e.m.substrates_count; i++)
	{
		dest.fraction_transferred_when_ingested()[i] = fraction_transferred_when_ingested()[i];
	}
}

void interactions_t::copy(interactions_t& dest)
{
	dest.apoptotic_phagocytosis_rate() = apoptotic_phagocytosis_rate();
	dest.necrotic_phagocytosis_rate() = necrotic_phagocytosis_rate();
	dest.other_dead_phagocytosis_rate() = other_dead_phagocytosis_rate();
	for (index_t i = 0; i < data_.e.cell_definitions_count; i++)
	{
		dest.live_phagocytosis_rates()[i] = live_phagocytosis_rates()[i];
	}
	for (index_t i = 0; i < data_.e.cell_definitions_count; i++)
	{
		dest.attack_rates()[i] = attack_rates()[i];
	}
	dest.attack_damage_rate() = attack_damage_rate();
	dest.attack_duration() = attack_duration();
	dest.total_damage_delivered() = total_damage_delivered();
	dest.attack_target() = attack_target();
	for (index_t i = 0; i < data_.e.cell_definitions_count; i++)
	{
		dest.immunogenicities()[i] = immunogenicities()[i];
	}
	for (index_t i = 0; i < data_.e.cell_definitions_count; i++)
	{
		dest.fusion_rates()[i] = fusion_rates()[i];
	}
}

void secretion_t::copy(secretion_t& dest)
{
	for (index_t i = 0; i < data_.e.m.substrates_count; i++)
	{
		dest.secretion_rates()[i] = secretion_rates()[i];
	}
	for (index_t i = 0; i < data_.e.m.substrates_count; i++)
	{
		dest.uptake_rates()[i] = uptake_rates()[i];
	}
	for (index_t i = 0; i < data_.e.m.substrates_count; i++)
	{
		dest.saturation_densities()[i] = saturation_densities()[i];
	}
	for (index_t i = 0; i < data_.e.m.substrates_count; i++)
	{
		dest.net_export_rates()[i] = net_export_rates()[i];
	}
}

void transformations_t::copy(transformations_t& dest)
{
	for (index_t i = 0; i < data_.e.cell_definitions_count; i++)
	{
		dest.transformation_rates()[i] = transformation_rates()[i];
	}
}

void volume_t::multiply_by_factor(real_t factor)
{
	total() *= factor;
	solid() *= factor;
	fluid() *= factor;

	nuclear() *= factor;
	nuclear_fluid() *= factor;
	nuclear_solid() *= factor;

	cytoplasmic() *= factor;
	cytoplasmic_fluid() *= factor;
	cytoplasmic_solid() *= factor;

	rupture_volume() *= factor;

	target_solid_cytoplasmic() *= factor;
	target_solid_nuclear() *= factor;
}

void volume_t::divide() { multiply_by_factor(0.5); }

void molecular_t::divide()
{
	for (index_t i = 0; i < data_.e.m.substrates_count; i++)
	{
		internalized_total_substrates()[i] *= 0.5;
	}
}

integrity_t::integrity_t(cell_data& data, const index_t& index) : phenotype_data_storage(data, index) {}

real_t& integrity_t::damage() { return data_.integrity.damages[index_]; }

real_t& integrity_t::damage_rate() { return data_.integrity.damage_rates[index_]; }

real_t& integrity_t::damage_repair_rate() { return data_.integrity.damage_repair_rates[index_]; }

void integrity_t::advance_damage(real_t time_step)
{
	damage() = (damage() + time_step * damage_rate()) / (1.0 + time_step * damage_repair_rate());
}

void integrity_t::set_defaults()
{
	damage() = 0.0;
	damage_rate() = 0.0;
	damage_repair_rate() = 0.0;
}

void integrity_t::copy(integrity_t& dest)
{
	dest.damage() = damage();
	dest.damage_rate() = damage_rate();
	dest.damage_repair_rate() = damage_repair_rate();
}
