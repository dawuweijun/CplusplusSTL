#include <opm/lattice/TwoPhaseLatticeBoltzmannSimulator.hpp>
#include <opm/lattice/GridManager.hpp>
#include <opm/lattice/FluidProperties.hpp>
#include <opm/lattice/LatticeBoltzmannModule.hpp>
#include <opm/lattice/LatticeBoltzmannSolver.hpp>
#include <opm/lattice/LatticeBoltzmannSolverOutput.hpp>
#include <opm/lattice/utility/StopWatch.hpp>
#include <opm/lattice/utility/SimulatorState.hpp>
#include <opm/lattice/utility/SimulatorTimer.hpp>
#include <opm/lattice/utility/updateState.hpp>
#include <boost/filesystem.hpp>

#include <numeric>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

class TwoPhaseLatticeBoltzmannSimulator::Impl
{
public:
    Impl(const GridManager& grid,
         const FluidProperties& red,
         const FluidProperties& blue,
         const LatticeBoltzmannModule& module);

    void
    run(SimulatorTimer& timer,
        SimulatorState& state);
private:
    bool output_;
    bool output_vtk_;
    std::string output_dir_;
    int output_interval_;
    const GridManager& grid_;
    const FluidProperties& red_;
    const FluidProperties& blue_;
    const LatticeBoltzmannModule& module_;
    //Solver.
    LatticeBoltzmannSolver solver_;
};


TwoPhaseLatticeBoltzmannSimulator::
TwoPhaseLatticeBoltzmannSimulator(const GridManager& grid,
                                  const FluidProperties& red,
                                  const FluidProperties& blue,
                                  const LatticeBoltzmannModule& module)
{
    pimpl_.reset(new Impl(grid, red, blue, module));
}

void
TwoPhaseLatticeBoltzmannSimulator::run(SimulatorTimer& timer, SimulatorState& state)
{
    return pimpl_->run(timer, state);
}

TwoPhaseLatticeBoltzmannSimulator::Impl::Impl(const GridManager& grid,
                                              const FluidProperties& red,
                                              const FluidProperties& blue,
                                              const LatticeBoltzmannModule& module)
    : grid_(grid)
    , red_(red)
    , blue_(blue)
    , module_(module)
    , solver_(grid, module, red, blue)
{
    output_ = true;
    if (output_) {
        output_vtk_ = true;
        output_dir_ = std::string("output");
        boost::filesystem::path fpath(output_dir_);
        try {
            create_directories(fpath);
        }
        catch (...) {
            std::cout << "Creating directories failed: " << fpath << std::endl;
            exit(1);
        }
        output_interval_ = 100;
    }
}

void
TwoPhaseLatticeBoltzmannSimulator::Impl::
run(SimulatorTimer& timer, SimulatorState& state)
{
    //main loop.
    StopWatch solver_timer;
    double stime = 0.0;
    StopWatch step_timer; 
    StopWatch total_timer;
    total_timer.start();
    std::fstream tstep_os;
    if (output_) {
        std::string filename = output_dir_ + "/step_timing.txt";
        tstep_os.open(filename.c_str(), std::fstream::out | std::fstream::app);
    }
    if (output_) {
        updateState(state, grid_, red_, blue_, module_);
        if (output_vtk_) {
            outputStateVtk(grid_, state, 0, output_dir_);
        }
        outputStateMatlab(grid_, state, 0, output_dir_);
    }
    while (!timer.done()) {
        step_timer.start();
        timer.report(std::cout);
        solver_timer.start();
        solver_.step(timer.currentStepLength(), state);
        solver_timer.stop();
        if (output_ && timer.currentStepNum()%output_interval_ == 0) {
            updateState(state, grid_, red_, blue_, module_);
            if (output_vtk_) {
                outputStateVtk(grid_, state, timer.currentStepNum()+1, output_dir_);
            }
            outputStateMatlab(grid_, state, timer.currentStepNum()+1, output_dir_);
        }
        const double st = solver_timer.secsSinceStart();
        std::cout << "\n     Lattice Boltzmann Simulator took: " << st << " seconds." << std::endl;
        stime += st;
        ++timer;
    }
    total_timer.stop();
    std::cout << "\n   Performance Summary:\n";
    std::cout << "\n     Total  Time taken: " << total_timer.secsSinceStart() << "s"
              << "\n     Solver Time taken: " << stime << "s";
}
