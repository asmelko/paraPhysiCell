#include <algorithm>
#include <array>
#include <iostream>

#include <BioFVM/types.h>

#include "simulator_measurer.h"

using namespace biofvm;
using namespace physicell;

#define make_duration_pair(D) std::make_pair(D, #D)

void simulator_durations::print_durations()
{
	std::array<std::pair<std::size_t, const char*>, 22> durations { make_duration_pair(diffusion),
																	make_duration_pair(secretion),
																	make_duration_pair(gradient),
																	make_duration_pair(host_sync),
																	make_duration_pair(device_sync),
																	make_duration_pair(custom_rules),
																	make_duration_pair(custom_interactions),
																	make_duration_pair(forces),
																	make_duration_pair(motility),
																	make_duration_pair(membrane),
																	make_duration_pair(spring),
																	make_duration_pair(position),
																	make_duration_pair(cell_cell),
																	make_duration_pair(container_mech),
																	make_duration_pair(mesh_mech),
																	make_duration_pair(neighbors_mech),
																	make_duration_pair(advance_phe),
																	make_duration_pair(container_phe),
																	make_duration_pair(mesh_phe),
																	make_duration_pair(neighbors_phe),
																	make_duration_pair(svg_save),
																	make_duration_pair(full_save) };

	std::size_t diff_total = 0;
	std::size_t mech_total = 0;
	std::size_t phen_total = 0;
	std::size_t save_total = 0;
	std::size_t total = 0;

	{
		index_t i = 0;
		for (; i < 4; i++)
			diff_total += durations[i].first;

		for (; i < 16; i++)
			mech_total += durations[i].first;

		for (; i < 20; i++)
			phen_total += durations[i].first;

		for (; i < 22; i++)
			save_total += durations[i].first;

		total = diff_total + mech_total + phen_total + save_total;
	}

	std::cout << "Duration ratios (total: " << (double)total / 1000 << "ms): "
			  << "Diffusion: " << (int)(diff_total / (double)total * 100)
			  << "% Mechanics: " << (int)(mech_total / (double)total * 100)
			  << "% Phenotype: " << (int)(phen_total / (double)total * 100)
			  << "% Save: " << (int)(save_total / (double)total * 100) << "%" << std::endl;

	std::sort(durations.begin(), durations.end(), [](auto& a, auto& b) { return a.first > b.first; });

	std::cout << "  ";
	for (index_t i = 0; i < 5; i++)
	{
		std::cout << durations[i].second << ": " << durations[i].first / 1000. << "ms ";
	}
	std::cout << std::endl;

	// zero them out
	diffusion = 0;
	secretion = 0;
	gradient = 0;
	custom_rules = 0;
	custom_interactions = 0;
	forces = 0;
	motility = 0;
	membrane = 0;
	spring = 0;
	position = 0;
	cell_cell = 0;
	container_mech = 0;
	mesh_mech = 0;
	neighbors_mech = 0;
	advance_phe = 0;
	container_phe = 0;
	mesh_phe = 0;
	neighbors_phe = 0;
	svg_save = 0;
	full_save = 0;
}
