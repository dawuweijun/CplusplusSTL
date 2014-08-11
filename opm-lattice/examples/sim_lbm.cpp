#include <opm/lattice/FluidProperties.hpp>
#include <opm/lattice/GridManager.hpp>
#include <opm/lattice/LatticeBoltzmannModule.hpp>
#include <opm/lattice/LatticeBoltzmannSolver.hpp>
#include <opm/lattice/TwoPhaseLatticeBoltzmannSimulator.hpp>
#include <opm/lattice/SimulatorState.hpp>
#include <opm/lattice/SimulatorTimer.hpp>
#include <opm/lattice/writeVtkData.hpp>
#include <opm/lattice/DataMap.hpp>
#include <opm/lattice/StopWatch.hpp>
#include <iostream>
#include <string>
#include <memory>

int main()
{
    std::cout << "\n================    Test program for two-phase lattice boltzmann simulator     ===============\n\n";
    std::shared_ptr<GridManager> grid;
    std::shared_ptr<FluidProperties> red;
    std::shared_ptr<FluidProperties> blue;
    std::shared_ptr<LatticeBoltzmannModule> module;
    SimulatorState state;
    SimulatorTimer simtimer;
    simtimer.init(1e6, 1);

    grid.reset(new GridManager(80, 20, 20));
    module.reset(new LatticeBoltzmannModule());
    red.reset(new FluidProperties(*grid, *module, 1.0, 1.0, 0.1667, 0.1));
    blue.reset(new FluidProperties(*grid, *module, 1.0, 1.0, 0.1667, 0.1));

    TwoPhaseLatticeBoltzmannSimulator simulator(*grid, *red, *blue, *module);
    std::cout << "\n\n================    Starting main simulation loop     ===============\n"
              << std::flush;

    simulator.run(simtimer, state);
    std::cout << "\n\n================    End of simulation     ===============\n\n";
    
    return 0;
}
