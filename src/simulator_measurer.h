#include <chrono>

namespace physicell {

struct simulator_durations
{
	std::chrono::high_resolution_clock::time_point start, end;

	std::size_t diffusion = 0;
	std::size_t secretion = 0;
	std::size_t gradient = 0;
	std::size_t host_sync = 0;
	std::size_t device_sync = 0;
	std::size_t custom_rules = 0;
	std::size_t custom_interactions = 0;
	std::size_t forces = 0;
	std::size_t motility = 0;
	std::size_t membrane = 0;
	std::size_t spring = 0;
	std::size_t position = 0;
	std::size_t cell_cell = 0;
	std::size_t container_mech = 0;
	std::size_t mesh_mech = 0;
	std::size_t neighbors_mech = 0;
	std::size_t advance_phe = 0;
	std::size_t container_phe = 0;
	std::size_t mesh_phe = 0;
	std::size_t neighbors_phe = 0;
	std::size_t full_save = 0;
	std::size_t svg_save = 0;

	void print_durations();
};

} // namespace physicell
