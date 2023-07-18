#include "interactions_solver.h"

#include <cmath>

#include "../../random.h"

using namespace biofvm;
using namespace physicell;

void ingest_volume(index_t lhs, index_t rhs, real_t* __restrict__ total, real_t* __restrict__ fluid,
				   real_t* __restrict__ solid, real_t* __restrict__ cytoplasmic, real_t* __restrict__ cytoplasmic_fluid,
				   real_t* __restrict__ cytoplasmic_solid, real_t* __restrict__ fluid_fraction,
				   real_t* __restrict__ cytoplasmic_to_nuclear_ratio, const real_t* __restrict__ nuclear_solid)
{
	// absorb fluid volume (all into the cytoplasm)
	cytoplasmic_fluid[lhs] += fluid[rhs];
	// absorb nuclear and cyto solid volume (into the cytoplasm)
	cytoplasmic_solid[lhs] += solid[rhs];

	// consistency calculations
	fluid[lhs] += fluid[rhs];
	solid[lhs] += solid[rhs];

	cytoplasmic[lhs] += fluid[rhs] + solid[rhs];
	total[lhs] += fluid[rhs] + solid[rhs];

	fluid_fraction[lhs] = fluid[lhs] / (total[lhs] + 1e-16);
	cytoplasmic_to_nuclear_ratio[lhs] = cytoplasmic_solid[lhs] / (nuclear_solid[lhs] + 1e-16);
}

void fuse_volume(index_t lhs, index_t rhs, real_t* __restrict__ total, real_t* __restrict__ fluid,
				 real_t* __restrict__ solid, real_t* __restrict__ cytoplasmic, real_t* __restrict__ cytoplasmic_fluid,
				 real_t* __restrict__ cytoplasmic_solid, real_t* __restrict__ nuclear,
				 real_t* __restrict__ nuclear_fluid, real_t* __restrict__ nuclear_solid,
				 real_t* __restrict__ fluid_fraction, real_t* __restrict__ cytoplasmic_to_nuclear_ratio,
				 real_t* __restrict__ target_solid_cytoplasmic, real_t* __restrict__ target_solid_nuclear)
{
	cytoplasmic_fluid[lhs] += cytoplasmic_fluid[rhs];
	nuclear_fluid[lhs] += nuclear_fluid[rhs];

	cytoplasmic_solid[lhs] += cytoplasmic_solid[rhs];
	nuclear_solid[lhs] += nuclear_solid[rhs];

	target_solid_cytoplasmic[lhs] += target_solid_cytoplasmic[rhs];
	target_solid_nuclear[lhs] += target_solid_nuclear[rhs];

	// consistency calculations
	fluid[lhs] = nuclear_fluid[lhs] + cytoplasmic_fluid[lhs];
	solid[lhs] = nuclear_solid[lhs] + cytoplasmic_solid[lhs];

	nuclear[lhs] = nuclear_fluid[lhs] + nuclear_solid[lhs];
	cytoplasmic[lhs] = cytoplasmic_fluid[lhs] + cytoplasmic_solid[lhs];
	total[lhs] = nuclear[lhs] + cytoplasmic[lhs];

	fluid_fraction[lhs] = fluid[lhs] / (total[lhs] + 1e-16);
	cytoplasmic_to_nuclear_ratio[lhs] = cytoplasmic_solid[lhs] / (nuclear_solid[lhs] + 1e-16);
}

void ingest_internalized(index_t lhs, index_t rhs, index_t substrates_count,
						 real_t* __restrict__ internalized_substrates,
						 const real_t* __restrict__ fraction_transferred_when_ingested)
{
	for (index_t i = 0; i < substrates_count; i++)
	{
		internalized_substrates[lhs * substrates_count + i] +=
			internalized_substrates[rhs * substrates_count + i]
			* fraction_transferred_when_ingested[rhs * substrates_count + i];
	}
}

void fuse_internalized(index_t lhs, index_t rhs, index_t substrates_count, real_t* __restrict__ internalized_substrates)
{
	for (index_t i = 0; i < substrates_count; i++)
	{
		internalized_substrates[lhs * substrates_count + i] += internalized_substrates[rhs * substrates_count + i];
	}
}

void fuse_position(index_t lhs, index_t rhs, index_t dims, real_t* __restrict__ total_volume,
				   real_t* __restrict__ position)
{
	const real_t both_volume = total_volume[lhs] + total_volume[rhs];

	for (int i = 0; i < dims; i++)
	{
		position[lhs * dims + i] =
			(position[lhs * dims + i] * total_volume[lhs] + position[rhs * dims + i] * total_volume[rhs]) / both_volume;
	}
}

void update_geometry(index_t i, real_t* __restrict__ radius, real_t* __restrict__ nuclear_radius,
					 real_t* __restrict__ surface_area, const real_t* __restrict__ total_volume,
					 const real_t* __restrict__ nuclear_volume)
{
	radius[i] = std::cbrt(total_volume[i] / (M_PI * 4.0 / 3.0));
	nuclear_radius[i] = std::cbrt(nuclear_volume[i] / (M_PI * 4.0 / 3.0));
	surface_area[i] = 4.0 * M_PI * radius[i] * radius[i];
}

void ingest(index_t lhs, index_t rhs, cell_data& data)
{
	ingest_volume(lhs, rhs, data.agent_data.volumes.data(), data.volumes.fluid.data(), data.volumes.solid.data(),
				  data.volumes.cytoplasmic.data(), data.volumes.cytoplasmic_fluid.data(),
				  data.volumes.cytoplasmic_solid.data(), data.volumes.fluid_fraction.data(),
				  data.volumes.cytoplasmic_to_nuclear_ratio.data(), data.volumes.nuclear_solid.data());

	ingest_internalized(lhs, rhs, data.e.m.substrates_count, data.agent_data.internalized_substrates.data(),
						data.agent_data.fraction_transferred_when_ingested.data());

	update_geometry(lhs, data.geometries.radius.data(), data.geometries.nuclear_radius.data(),
					data.geometries.surface_area.data(), data.agent_data.volumes.data(), data.volumes.nuclear.data());

	data.to_remove[rhs] = 1;
}

void fuse(index_t lhs, index_t rhs, cell_data& data)
{
	fuse_position(lhs, rhs, data.e.mechanics_mesh.dims, data.agent_data.volumes.data(),
				  data.agent_data.positions.data());

	fuse_volume(lhs, rhs, data.agent_data.volumes.data(), data.volumes.fluid.data(), data.volumes.solid.data(),
				data.volumes.cytoplasmic.data(), data.volumes.cytoplasmic_fluid.data(),
				data.volumes.cytoplasmic_solid.data(), data.volumes.nuclear.data(), data.volumes.nuclear_fluid.data(),
				data.volumes.nuclear_solid.data(), data.volumes.fluid_fraction.data(),
				data.volumes.cytoplasmic_to_nuclear_ratio.data(), data.volumes.target_solid_cytoplasmic.data(),
				data.volumes.target_solid_nuclear.data());

	fuse_internalized(lhs, rhs, data.e.m.substrates_count, data.agent_data.internalized_substrates.data());

	update_geometry(lhs, data.geometries.radius.data(), data.geometries.nuclear_radius.data(),
					data.geometries.surface_area.data(), data.agent_data.volumes.data(), data.volumes.nuclear.data());

	data.number_of_nuclei[lhs] += data.number_of_nuclei[rhs];

	data.to_remove[rhs] = 1;
}

void attack(index_t lhs, index_t rhs, real_t time_step, const real_t* __restrict__ damage_rate,
			real_t* __restrict__ damage, real_t* __restrict__ total_attack_time)
{
	damage[rhs] += damage_rate[lhs] * time_step;
	total_attack_time[rhs] += time_step;
}

void update_cell_cell_interactions_internal(
	index_t n, index_t cell_defs_count, real_t time_step, const std::uint8_t* __restrict__ dead,
	const index_t* __restrict__ cell_definition_indices, const real_t* __restrict__ dead_phagocytosis_rate,
	const real_t* __restrict__ live_phagocytosis_rate, const real_t* __restrict__ attack_rate,
	const real_t* __restrict__ fusion_rate, const real_t* __restrict__ immunogenicity,
	const std::vector<index_t>* __restrict__ neighbors, const std::uint8_t* __restrict__ to_remove, cell_data& data)
{
	for (index_t cell_index = 0; cell_index < n; cell_index++)
	{
		bool attacked_once = false;
		bool phagocytosed_once = false;
		bool fused_once = false;

		if (dead[cell_index] == 1 || to_remove[cell_index] == 1)
		{
			continue;
		}

		for (std::size_t i = 0; i < neighbors[cell_index].size(); i++)
		{
			const index_t neighbor_index = neighbors[cell_index][i];

			if (to_remove[neighbor_index] == 1)
			{
				continue;
			}

			const index_t neighbor_cell_def_index = cell_definition_indices[neighbor_index];

			const auto random_number = random::instance().uniform();

			if (dead[neighbor_index] == 1)
			{
				if (random_number < dead_phagocytosis_rate[cell_index] * time_step)
				{
					ingest(cell_index, neighbor_index, data);
				}

				continue;
			}

			if (!phagocytosed_once
				&& random_number
					   < live_phagocytosis_rate[cell_index * cell_defs_count + neighbor_cell_def_index] * time_step)
			{
				ingest(cell_index, neighbor_index, data);

				phagocytosed_once = true;

				continue;
			}

			{
				const auto attack_r = attack_rate[cell_index * cell_defs_count + neighbor_cell_def_index];
				const auto immuno_r =
					immunogenicity[neighbor_index * cell_defs_count + cell_definition_indices[cell_index]];

				if (!attacked_once && random_number < attack_r * immuno_r * time_step)
				{
					attack(cell_index, neighbor_index, time_step, data.damage.data(),
						   data.interactions.damage_rate.data(), data.total_attack_time.data());

					attacked_once = true;

					continue;
				}
			}

			if (!fused_once
				&& random_number < fusion_rate[cell_index * cell_defs_count + neighbor_cell_def_index] * time_step)
			{
				fuse(cell_index, neighbor_index, data);

				fused_once = true;

				continue;
			}
		}
	}
}

void interactions_solver::update_cell_cell_interactions(environment& e)
{
	auto& data = get_cell_data(e);

	update_cell_cell_interactions_internal(
		data.agents_count, e.cell_definitions_count, e.mechanics_time_step, data.death.dead.data(),
		data.cell_definition_indices.data(), data.interactions.dead_phagocytosis_rate.data(),
		data.interactions.live_phagocytosis_rates.data(), data.interactions.attack_rates.data(),
		data.interactions.fussion_rates.data(), data.interactions.immunogenicities.data(), data.neighbors.data(),
		data.to_remove.data(), data);
}
