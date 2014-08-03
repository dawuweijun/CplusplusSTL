#include <opm/lattice/FluidProperties.hpp>
#include <opm/lattice/GridManager.hpp>
#include <opm/lattice/LatticeBoltzmannModule.hpp>
#include <opm/lattice/LatticeBoltzmannSolver.hpp>
#include <opm/lattice/TwoPhaseLatticeBoltzmannSimulator.hpp>
#include <opm/lattice/SimulatorState.hpp>
#include <opm/lattice/writeVtkData.hpp>
#include <opm/lattice/SimulatorTimer.hpp>
#include <opm/lattice/DataMap.hpp>
#include <boost/filesystem.hpp>
#include <opm/lattice/StopWatch.hpp>
#include <iostream>
#include <string>

int main()
{
    std::cout << "\n================    Test program for two-phase lattice boltzmann simulator     ===============\n\n";
    
    
    const int NX = 80;
    const int NY = 20;
    const int NZ = 20;
    std::shared_ptr<GridManager> grid;
    std::shared_ptr<FluidProperties> red;
    std::shared_ptr<FluidProperties> blue;
    std::shared_ptr<LatticeBoltzmannModule> module;
    grid.reset(new GridManager(NX, NY, NZ));
    module.reset(new LatticeBoltzmannModule());
    red.reset(new FluidProperties(*grid, *module, 1.0, 1.0, 0.1667, 0.1));
    blue.reset(new FluidProperties(*grid, *module, 1.0, 1.0, 0.1667, 0.1));
    /*
    GridManager grid(80, 20, 20);
    LatticeBoltzmannModule module();
    FluidProperties red(grid, module, 1.0, 1.0, 0.1667, 0.1);
    FluidProperties blue(grid, module, 1.0, 1.0, 0.1667, 0.1);
    */
    SimulatorState state;


    std::cout << "\n\n================    Starting main simulation loop     ===============\n"
              << std::flush;
    SimulatorTimer simtimer;

    TwoPhaseLatticeBoltzmannSimulator simulator(*grid, *red, *blue, *module);
    simtimer.init(1e6, 1);
    simulator.run(simtimer, state);
    std::cout << "\n\n================    End of simulation     ===============\n\n";
    
    return 0;
}
