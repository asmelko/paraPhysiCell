#include "cell.h"

#include <BioFVM/microenvironment.h>

#include "cell_container.h"
#include "cell_data.h"

using namespace biofvm;
using namespace physicell;

cell::cell(agent_id_t id, cell_data& data, index_t index) : agent(id, data.agent_data, index), data_(data) {}

real_t* cell::velocities() { return data_.velocities.data() + index_ * data_.m.mesh.dims; }
