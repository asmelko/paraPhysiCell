# paraPhysiCell: the optimized parallel implementation of PhysiCell core features

This repository contains the parallel CPU implementation of PhysiCell. The speedups vary from 2x to 5x for small 2D simulations and ~60x on big 3D simulations. The BioFVM part is implemented using the optimized [paraBioFVM library](https://github.com/asmelko/paraBioFVM). The goal of this project is twofold:
 - implement optimized cell mechanics (their interactions, forces, motility, spring attachments, ...)
 - and propose the architectural changes for better sustainability and extendability of the code.

> The project capabilities allow paraPhysiCell to be used as a reasonable substitute for the original PhysiCell in *PhysiCell Studio*.

The missing features (not a complete list): 
- anything intracellular (PhysiBoSS, libRoadrunner, dFBA),
- crazy make targets, such as make-movie or make-zip, 
- PhysiMeSS.

The version of PhysiCell this project builds upon: **1.12.0**.

## Build requirements

The only requirements are CMake and C++20 compliant compiler. To build everything (all the libraries, unit tests, sample projects, ...), write:
```bash
# configure cmake for release build in 'build' directory:
cmake -DCMAKE_BUILD_TYPE=Release -B build 
# build the contents of 'build' directory:
cmake --build build 
```
The sample projects can be then accessed in `build/sample_projects/<project>` directory.

For build in unix-like environments, `gcc` and `clang` were tested and both are supported (the preferred one is gcc, since it seems to generate much more optimized code). For build under Windows, only `clang-cl` compiler is currently supported (apart from MinGW or WSL options).

Optional: For the ability to run in parallel, the compiler with OpenMP support is required. 

## Optimizations

The optimizations in paraPhysiCell are solely in mechanics timestep. It encouples updating cell velocities, running cell-cell interactions and handling operations on the cell container. 

To achieve more drastic optimizations, some restrictions have been made on paraPhysiCell. These are mainly:
- the `update_velocity`, `calculate_distance_to_membrane` and `add_cell_basement_membrane_interactions` cell functions are now globally set to be same for each cell,
- `use_2D` flag in motility phenotype is now ignored. We are either simulating in 2D and `use_2D` is globally true or we are simulating in 3D and `use_2D` is globally false,
- now, when simulating 2D simulation, the simulation is truly 2D --- meaning that cell position stores only 2 coordinates. 

### GPU support
The project makes use of paraBioFVM's GPU support. It can be simply enabled by passing `-DBUILD_FOR_DEVICE` flag to cmake configure step. However, it may be worth it only for large 3D simulations due to its current suboptimalities. For more info see [paraBioFVM README](https://github.com/asmelko/paraBioFVM).

## Architectural changes

There has been some changes in the way how simulation code is written, how the global environment is configured and how the directory structures is composed in the pursuit to simplify the life of a PhysiCell user.

### Directory structure

Since this is just a light version of PhysiCell, the directory structure is very simple. In the top directory, we can find:
- `sample_projects` which contains sample projects, 
- `src` which contains all the source files needed to compile the "core" of the paraPhysiCell (the CMake target is therefore not-so-surprisingly named *PhysiCellCore*), 
- `test` which contains unit tests.

`src` divides the ordinary sources from the most performant critical ones, which are stored in `solver` directory. Lastly, `src/original` directory contains **the original PhysiCell sources** which were required for this project, namely the ones from the original `core` and `modules` directories. 
> The sources there are formated using `clang-format` and have their data types changed to paraPhysiCell ones. Apart from that, the sources are almost untouched. **Hence, in case of new PhysiCell release, these files must be checked and updated appropriately.** 

### New classes

This project was implemented with the aim to be very flexible with respect to the future extensions. To achieve so, we targeted to improve in these 3 domains:

#### Encapsule the whole PhysiCell state into one class. 

This is important to simplify the reasoning about the code structure, removing the hidden dependencies and defining the single source of truth. For that, we implemented the structure `environment`, which is simply the union of BioFVM's microenvironment, PhysiCell's mechanics mesh, BioFVM's agent container and PhysiCell's cell definitions (rules and signals are still missing though). 

Same as in paraBioFVM's `microenvironment_builder`, we implemented the class `builder`, which contains convenience functions for customizing the environment object. 

#### Make the code extensible wrt new simulation steps.

The original sample projects contain a big amount of repeating boilerplate code just to make the simulation run. We alleviate this burden by implementing the `simulator` class, which contains all the simulation code for diffusion/mechanics/phenotype timesteps and provides the interface for specifying the custom ones. The simplest example of a custom timestep would look like this:
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

#### Make the code extensible wrt new cell types.

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

If one wants the default type of container to be something different than `cell` (for example to have cells with legs as a default), they can simply assign a new specialization of `agent_container_common` to environment cell container by
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
Go have a look at the some reimplemented sample projects!
