# paraPhysiCell: the optimized parallel implementation of PhysiCell core features

```
general:
- elements count (N) = detachment_rate
- cell radius (R_cell) = specified in custom.cpp
- viscosity (n) = attachment_rate
- equilibrium distance (r_eq) = cell_cell_repulsion_strength (automatically computed in custom.cpp)
- element_radius (R_elem) = radius (automatically computed in custom.cpp)



intra:
- scaling_factor (rho) = relative_maximum_adhesion_distance
- stiffness (k) = attachment_elastic_constant

inter:
- scaling_factors (rho) = live_phagocytosis_rate
- stiffness (k) = attack_rates
- equilibrium multiplier = cell_adhesion_affinity (when =1, it is just sum of 2 element radii)

1. build the project as is specified below
 - the binaries and configs will be in ./build/sample_projects/(sem-test/sem-playground)
2. run physicell studio on sem-test or sem-playground

```

This repository contains the parallel CPU implementation of PhysiCell. The BioFVM part is implemented using the optimized [paraBioFVM library](https://github.com/asmelko/paraBioFVM). The goal of this project is twofold:
 - implement optimized cell mechanics (their interactions, forces, motility, spring attachments, ...)
 - and propose architectural changes for better sustainability and extendability of the code.

> The project capabilities allow paraPhysiCell to be used as a reasonable substitute for the original PhysiCell in *PhysiCell Studio*.

The missing features (not a complete list):
- anything intracellular (PhysiBoSS, libRoadrunner, dFBA),
- crazy make targets, such as make-movie or make-zip,
- PhysiMeSS.

The version of PhysiCell this project builds upon: **1.12.0**.

## Simple Benchmarks

This table shows the wall times of a modified `heterogeneity-sample` simulation. The first row shows the simulation modifications and the first column shows the hardware the simulation ran on. The first wall time is for the optimized paraPhysiCell and the second wall time in the parentheses is the original Physicell. **These are very simple and naive measurements. They are here just to show the achievable speedups. Please, take them with a grain of salt.**

 | | 100x100 voxels, 900 agents, 600m simul. time | 300x300 voxels, 14300 agents, 600m simul. time | 100x100x100 voxels, 14300 agents, 60m simul. time | 200x200x200 voxels, 88700 agents, 60m simul. time |
 | - | - | - | - | - |
 | i7-8650U laptop CPU (4 HT cores) | 4.5s (12.7s) | 38.3s (256.3s) | 18.6s (246.3s) | - |
 | 2 sockets Intel Xeon Platinum 8160 CPU (2x24 cores) | 18.7s (14.4s) | 27.4s (116.9s) | 6.8s (156.6s) | 36.8s (1231.6s) |


## Build requirements

The only requirements are git, CMake and C++20 compliant compiler. To build everything (all the libraries, unit tests, sample projects, ...), write:
```bash
# Fetch dependencies
git submodule update --init --recursive
# configure cmake for release build in 'build' directory:
cmake -DCMAKE_BUILD_TYPE=Release -B build
# build the contents of 'build' directory:
cmake --build build
```
The sample projects can be then accessed in `build/sample_projects/<project>` directory.

For a build in Unix-like environments, `gcc` and `clang` were tested and both are supported (the preferred one is gcc, since it seems to generate much more optimized code). For a build under Windows, only `clang-cl` compiler is currently supported (apart from MinGW or WSL options).

Optional: For the ability to run in parallel, the compiler with OpenMP support is required.

## Optimizations

The optimizations in paraPhysiCell are solely in mechanics timestep. It encouples updating cell velocities, running cell-cell interactions and handling operations on the cell container.

To achieve more drastic optimizations, some restrictions have been made on paraPhysiCell. These are mainly:
- the `update_velocity`, `calculate_distance_to_membrane` and `add_cell_basement_membrane_interactions` cell functions are now globally set to be the same for each cell,
- `use_2D` flag in motility phenotype is now ignored. We are either simulating in 2D and `use_2D` is globally true or we are simulating in 3D and `use_2D` is globally false,
- now, when simulating 2D simulation, the simulation is truly 2D --- meaning that cell position stores only 2 coordinates.

### GPU support
The project makes use of paraBioFVM's GPU support. It can be simply enabled by passing `-DBUILD_FOR_DEVICE` flag to cmake configure step. However, it may be worth it only for large 3D simulations due to its current suboptimalities. For more info see [paraBioFVM README](https://github.com/asmelko/paraBioFVM).

## Architectural changes

There have been some changes in the way simulation code is written, how the global environment is configured and how the directory structures are composed in the pursuit to simplify the life of a PhysiCell user.

### Directory structure

Since this is just a light version of PhysiCell, the directory structure is very simple. In the top directory, we can find:
- `sample_projects` which contains sample projects,
- `src` which contains all the source files needed to compile the "core" of the paraPhysiCell (the CMake target is therefore not-so-surprisingly named *PhysiCellCore*),
- `test` which contains unit tests.

`src` divides the ordinary sources from the most performant critical ones, which are stored in `solver` directory. Lastly, `src/original` directory contains **the original PhysiCell sources** which were required for this project, namely the ones from the original `core` and `modules` directories.
> The sources there are formated using `clang-format` and have their data types changed to paraPhysiCell ones. Apart from that, the sources are almost untouched. **Hence, in case of new PhysiCell release, these files must be checked and updated appropriately.**

### New classes

This project was implemented with the aim to be very flexible with respect to future extensions. To achieve so, we targeted to improve in these 3 domains:

#### 1. Encapsule the whole PhysiCell state into one class.

This is important to simplify the reasoning about the code structure, removing the hidden dependencies and defining the single source of truth. For that, we implemented the structure `environment`, which is simply the union of BioFVM's microenvironment, PhysiCell's mechanics mesh, BioFVM's agent container and PhysiCell's cell definitions (rules and signals are still missing though).

Same as in paraBioFVM's `microenvironment_builder`, we implemented the class `builder`, which contains convenience functions for customizing the environment object.

#### 2. Make the code extensible wrt. new simulation steps.

The original sample projects contain a large amount of repeating boilerplate code just to make the simulation run. We alleviate this burden by implementing the `simulator` class, which contains all the simulation code for diffusion/mechanics/phenotype timesteps and provides the interface for specifying the custom ones. The simplest example of a custom timestep would look like this:
```c++
simulation s;
real_t pension_time_step = 60'000;

// this adds the new code inside the simulation loop
// running once every 60'000 simulated minutes
s.add_simulation_step(pension_time_step, [](environment& e)
{
    // make random cell retire
});

// runs simulation with the standard diffusion/mechanics/phenotype steps + new pension step
s.run(...);
```

#### 3. Make the code extensible wrt. new cell types.

While the memory layout of the cell data has changed, the `cell` class and the `cell_container` are implemented in a way to allow a reasonable freedom of inheritance and polymorphism. Take a new cell class as an example:
```c++
class cell_with_legs : public cell
{
    leg_t left, right;

    cell_with_legs(...) {/* a custom constructor here */}
};
```
This cell can be instanciated simply by calling cell container's `create`
```c++
env.get_container().create<cell_with_legs>();
```
which would now add an object of type `cell_with_legs` inside the array of objects that are primarily of type `cell`.

If one wants the default type of container to be something different than `cell` (for example to have cells with legs as a default), they can simply assign a new specialization of `agent_container_common` to the environment cell container by
```c++
using better_container = biofvm::agent_container_common<cell_with_legs, cell_data>;
env.m.agents = std::make_unique<better_container>();
```
. Then, each call to `env.get_container<better_container>().create()` would create a `cell_with_legs` class by default.

---

Thanks to these architectural changes, one can write the custom sample project faster with doors very well open for extension. The simplest sample project code can now look like this:
```c++
int main(int argc, char* argv[])
{
    // create builder
    physicell::builder builder(argc, argv);

    /* configure env & microenv */

    // builds environment using the config file passed in args
    // i.e., builds microenv mesh, container, cell definitions, rules, signals
    auto e = builder.build_environment();

    // create simulator
    simulator s(*e);

    // initialize and run it with PhysiCell_Settings
    s.initialize(builder.get_settings());
    s.run(builder.get_settings(), some_coloring_func);
}
```
Go have a look at some reimplemented sample projects!
