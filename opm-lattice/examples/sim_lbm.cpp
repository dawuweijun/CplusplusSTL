#include <opm/lattice/FluidProperties.hpp>
#include <opm/lattice/GridManager.hpp>
#include <opm/lattice/LatticeBoltzmannModule.hpp>
#include <opm/lattice/LatticeBoltzmannSolver.hpp>
#include <opm/lattice/TwoPhaseLatticeBoltzmannSimulator.hpp>
#include <opm/lattice/utility/SimulatorState.hpp>
#include <opm/lattice/utility/SimulatorTimer.hpp>
#include <opm/lattice/utility/initState.hpp>
#include <iostream>
#include <ostream>
#include <iterator>
#include <string>
#include <memory>
#include <algorithm>

int main()
{
    std::cout << "\n================    Test program for two-phase lattice boltzmann simulator     ===============\n\n";
    std::shared_ptr<GridManager> grid;
    std::shared_ptr<FluidProperties> red;
    std::shared_ptr<FluidProperties> blue;
    std::shared_ptr<LatticeBoltzmannModule> module;
    SimulatorState state;
    SimulatorTimer simtimer;
    simtimer.init(7000, 1);

    grid.reset(new GridManager(80, 20, 20));
    module.reset(new LatticeBoltzmannModule());
    red.reset(new FluidProperties(*module, 1.0, 1.0, 0.1667, 0.1));
    blue.reset(new FluidProperties(*module, 1.0, 1.0, 0.1667, 0.1));
 //   blue.reset(new FluidProperties(*module, 0, 1., 0., 0.));

    red->setPotential(0.0, -0.281, -0.14, -0.01);
    blue->setPotential(0.0, -0.281, -0.14, 0.01);

    state.init(*grid, *module);
    initState(state, *grid, *red, *blue, *module);
 
    TwoPhaseLatticeBoltzmannSimulator simulator(*grid, *red, *blue, *module);

    std::cout << "\n\n================    Starting main simulation loop     ===============\n"
              << std::flush;
    simulator.run(simtimer, state);
    std::cout << "\n\n================    End of simulation     ===============\n\n";
    
    return 0;
}
