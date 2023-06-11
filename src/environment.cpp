#include "environment.h"

using namespace biofvm;
using namespace physicell;

environment::environment(microenvironment& m) : m(m) {}

cell_container_base& environment::cells() { return cast_container<cell_container_base&>(); }
