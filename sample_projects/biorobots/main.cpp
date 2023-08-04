#include <exception>

#include "custom.h"
#include "src/simulator.h"

using namespace biofvm;
using namespace physicell;

int main(int argc, char* argv[])
{
	try
	{
		physicell::builder builder(argc, argv);

		setup_microenvironment(builder.get_microenvironment_builder());

		create_cell_types(builder);

		auto e = builder.build_environment();

		setup_tissue(*e, builder.get_parameters(), builder.get_config_root());

		simulator s;

		s.initialize(*e);

		s.run(*e, builder.get_settings(), get_robot_coloring_function(builder.get_parameters()));
	}
	catch (std::exception& e)
	{
		std::cerr << "Exception caught: " << e.what() << std::endl;
		return 1;
	}

	return 0;
}
