#include "environment.h"

#include "models/standard_position_model.h"

using namespace biofvm;
using namespace physicell;

environment::environment(microenvironment& m, index_t cell_definitions_count,
						 biofvm::point_t<biofvm::index_t, 3> mechanics_voxel_shape)
	: m(m),
	  virtual_wall_at_domain_edges(false),
	  rules_enabled(false),
	  automated_spring_adhesion(true),
	  divisions_count(0),
	  deaths_count(0),
	  mechanics_mesh(m.mesh.dims, m.mesh.bounding_box_mins, m.mesh.bounding_box_maxs, mechanics_voxel_shape),
	  mechanics_time_step(0.1),
	  phenotype_time_step(6),
	  current_time(0),
	  cell_definitions_count(cell_definitions_count),
	  cell_definitions_data(*this),
	  position(std::make_unique<standard_position_model>())
{
	cells_in_mechanics_voxels = std::make_unique<std::vector<index_t>[]>(mechanics_mesh.voxel_count());
}

cell_container_base& environment::container_base() { return get_container<cell_container_base&>(); }

cell_definition& environment::create_cell_definition()
{
	cell_definitions_data.add();
	cell_definitions.emplace_back(std::make_unique<cell_definition>(*this, cell_definitions.size()));
	return *cell_definitions.back();
}

cell_definition& environment::cell_defaults() { return *cell_definitions[0]; }

cell_definition* environment::find_cell_definition(const std::string& name)
{
	auto it = std::find_if(cell_definitions.begin(), cell_definitions.end(),
						   [&name](const std::unique_ptr<cell_definition>& cd) { return cd->name == name; });

	return it == cell_definitions.end() ? nullptr : it->get();
}

cell_definition* environment::find_cell_definition(index_t type)
{
	auto it = std::find_if(cell_definitions.begin(), cell_definitions.end(),
						   [&type](const std::unique_ptr<cell_definition>& cd) { return cd->type == type; });

	return it == cell_definitions.end() ? nullptr : it->get();
}

void environment::display_info()
{
	m.display_info();
	display_cell_definitions_info();
}

void environment::display_cell_definitions_info()
{
	for (std::size_t n = 0; n < cell_definitions.size(); n++)
	{
		cell_definition* pCD = cell_definitions[n].get();
		std::cout << n << " :: type:" << pCD->type << " name: " << pCD->name << std::endl;

		// summarize cycle model
		if (pCD->phenotype.cycle.pCycle_Model != NULL)
		{
			std::cout << "\t cycle model: " << pCD->phenotype.cycle.model().name
					  << " (code=" << pCD->phenotype.cycle.model().code << ")" << std::endl;

			// let's show the transition rates
			Cycle_Model* pCM = &(pCD->phenotype.cycle.model());
			Cycle_Data* pCMD = &(pCD->phenotype.cycle.data);
			for (std::size_t n = 0; n < pCM->phases.size(); n++)
			{
				std::cout << "\t\tPhase " << n << ": " << pCM->phases[n].name << std::endl;
			}
			std::cout << "\t\tCycle transitions: " << std::endl << "\t\t-----------------" << std::endl;
			for (std::size_t n = 0; n < pCM->phase_links.size(); n++)
			{
				for (std::size_t k = 0; k < pCM->phase_links[n].size(); k++)
				{
					int start = pCM->phase_links[n][k].start_phase_index;
					int end = pCM->phase_links[n][k].end_phase_index;
					std::cout << "\t\t" << pCM->phases[start].name << " --> " << pCM->phases[end].name
							  << " w mean duration " << 1.0 / pCMD->transition_rate(start, end) << " min" << std::endl;
				}
			}
		}
		else
		{
			std::cout << "\t cycle model not initialized" << std::endl;
		}

		// summarize death models
		std::cout << "\t death models: " << std::endl;
		for (std::size_t k = 0; k < pCD->phenotype.death.models.size(); k++)
		{
			std::cout << "\t\t" << k << " : " << pCD->phenotype.death.models[k]->name
					  << " (code=" << pCD->phenotype.death.models[k]->code << ")"
					  << " with rate " << pCD->phenotype.death.rates[k] << " 1/min" << std::endl;

			Cycle_Model* pCM = (pCD->phenotype.death.models[k]);
			Cycle_Data* pCMD = &(pCD->phenotype.death.models[k]->data);


			std::cout << "\t\tdeath phase transitions: " << std::endl << "\t\t------------------------" << std::endl;
			for (std::size_t n = 0; n < pCM->phase_links.size(); n++)
			{
				for (std::size_t k = 0; k < pCM->phase_links[n].size(); k++)
				{
					int start = pCM->phase_links[n][k].start_phase_index;
					int end = pCM->phase_links[n][k].end_phase_index;
					std::cout << "\t\t" << pCM->phases[start].name << " --> " << pCM->phases[end].name
							  << " w mean duration " << 1.0 / pCMD->transition_rate(start, end) << " min" << std::endl;
				}
			}
		}

		auto bool_to_str = [](bool b) { return b ? "true" : "false"; };

		// summarize functions
		cell_functions* pCF = &(pCD->functions);
		std::cout << "\t key functions: " << std::endl;
		std::cout << "\t\t migration bias rule: ";
		std::cout << bool_to_str(pCD->phenotype.motility.update_migration_bias_direction() != nullptr);
		std::cout << std::endl;
		std::cout << "\t\t custom rule: ";
		std::cout << bool_to_str(pCF->custom_cell_rule != nullptr);
		std::cout << std::endl;
		std::cout << "\t\t phenotype rule: ";
		std::cout << bool_to_str(pCF->update_phenotype != nullptr);
		std::cout << std::endl;
		std::cout << "\t\t volume update function: ";
		std::cout << bool_to_str(pCF->volume_update_function != nullptr);
		std::cout << std::endl;
		std::cout << "\t\t mechanics function: ";
		std::cout << "true";
		std::cout << std::endl;
		std::cout << "\t\t contact function: ";
		std::cout << bool_to_str(pCF->contact_function != nullptr);
		std::cout << std::endl;

		// summarize motility

		motility_t* pM = &(pCD->phenotype.motility);
		std::string val = "true";
		if (pM->is_motile() == false)
		{
			val = "false";
		}

		std::string dimen = "3D";
		if (m.mesh.dims == 2)
		{
			dimen = "2D";
		}

		std::cout << "\tmotility (enabled: " << val << " in " << dimen << ")" << std::endl
				  << "\t\tspeed: " << pM->migration_speed() << " micron/min" << std::endl
				  << "\t\tbias: " << pM->migration_bias() << " " << std::endl
				  << "\t\tpersistence time: " << pM->persistence_time() << " min" << std::endl
				  << "\t\tchemotaxis (enabled: ";

		val = "maybe";
		if (pM->update_migration_bias_direction() == nullptr)
		{
			val = "no";
		}
		std::cout << val << ")" << std::endl
				  << "\t\t\talong " << pM->chemotaxis_direction() << " * grad("
				  << m.substrates_names[pM->chemotaxis_index()] << ") " << std::endl;

		// secretion



		// mechanics

		mechanics_t* pMech = &(pCD->phenotype.mechanics);

		std::cout << "\tmechanics:" << std::endl
				  << "\t\tcell_cell_adhesion_strength: " << pMech->cell_cell_adhesion_strength() << std::endl
				  << "\t\tcell_cell_repulsion_strength: " << pMech->cell_cell_repulsion_strength() << std::endl
				  << "\t\trel max adhesion dist: " << pMech->relative_maximum_adhesion_distance() << std::endl
				  << "\t\tcell_BM_adhesion_strength: " << pMech->cell_BM_adhesion_strength() << std::endl
				  << "\t\tcell_BM_repulsion_strength: " << pMech->cell_BM_repulsion_strength() << std::endl
				  << "\t\tattachment_elastic_constant: " << pMech->attachment_elastic_constant() << std::endl
				  << "\t\tattachment_rate: " << pMech->attachment_rate() << std::endl
				  << "\t\tdetachment_rate: " << pMech->detachment_rate() << std::endl;

		// size


		// intracellular
		// if (pCD->phenotype.intracellular != NULL)
		// {
		// 	pCD->phenotype.intracellular->display(std::cout);
		// }

		Custom_Cell_Data* pCCD = &(pCD->custom_data);
		std::cout << "\tcustom data: " << std::endl;
		for (std::size_t k = 0; k < pCCD->variables.size(); k++)
		{
			std::cout << "\t\t" << pCCD->variables[k] << std::endl;
		}
		std::cout << "\tcustom vector data: " << std::endl;
		for (std::size_t k = 0; k < pCCD->vector_variables.size(); k++)
		{
			std::cout << "\t\t" << pCCD->vector_variables[k] << std::endl;
		}
		std::cout << "\t\t\tNOTE: custom vector data will eventually be merged with custom data" << std::endl;
	}

	return;
}
