#include <opm/lattice/TwoPhaseLatticeBoltzmannSimulator.hpp>
#include <boost/filesystem.hpp>

#include <numeric>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>


TwoPhaseLatticeBoltzmannSimulator::
TwoPhaseLatticeBoltzmannSimulator(const GridManager& grid,
                                  const FluidProperties& red,
                                  const FluidProperties& blue,
                                  const LatticeBoltzmannModule& module)
    : grid_(grid)
    , red_(red)
    , blue_(blue)
    , module_(module)
    , solver_(grid_, module_, red_, blue_)
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
TwoPhaseLatticeBoltzmannSimulator::
outputStateVtk(const GridManager& grid,
               const SimulatorState& state,
               const int step,
               const std::string& output_dir)
{
    std::ostringstream vtkfilename;
    vtkfilename << output_dir << "/vtk_files";
    boost::filesystem::path fpath(vtkfilename.str());
    try {
        create_directories(fpath);
    }
    catch (...) {
        std::cout << "Creating directories failed: " << fpath << std::endl;
        exit(1);
    }
    vtkfilename << "/output-" << std::setw(6) << std::setfill('0') << step << ".vtu";
    std::ofstream vtkfile(vtkfilename.str().c_str());
    if (!vtkfile) {
        std::cout << "Failed to open " << vtkfilename.str() << std::endl;
        exit(1);
    }
    DataMap dm;
    dm["density"] = &state.velocity();
    dm["pressure"] = &state.pressure();
    writeVtkData(grid, dm, vtkfile);
}


void 
TwoPhaseLatticeBoltzmannSimulator::
outputStateMatlab(const GridManager& grid,
                  const SimulatorState& state,
                  const int step,
                  const std::string& output_dir)
{
    DataMap dm;
    dm["density"] = &state.velocity();
    dm["pressure"] = &state.pressure();
    // Write data (not grid) in Matlab format
    for (DataMap::const_iterator it = dm.begin(); it != dm.end(); ++it) {
        std::ostringstream fname;
        fname << output_dir << "/" << it->first;
        boost::filesystem::path fpath = fname.str();
        try {
            create_directories(fpath);
        }
        catch (...) {
            std::cout << "Creating directories failed: " << fpath << std::endl;
            exit(1);
        }
        fname << "/" << std::setw(6) << std::setfill('0') << step << ".txt";
        std::ofstream file(fname.str().c_str());
        if (!file) {
            std::cout << "Failed to open " << fname.str() << std::endl;
            exit(1);
        }
        file.precision(15);
        const std::vector<double>& d = *(it->second);
        std::copy(d.begin(), d.end(), std::ostream_iterator<double>(file, "\n"));
    }
}

void
TwoPhaseLatticeBoltzmannSimulator::
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
    while (!timer.done()) {
        step_timer.start();
        timer.report(std::cout);
        if (output_ && (timer.currentStepNum() % output_interval_ == 0)) {
            if (output_vtk_) {
                outputStateVtk(grid_, state, timer.currentStepNum(), output_dir_);
            }
            outputStateMatlab(grid_, state, timer.currentStepNum(), output_dir_);
        }
        solver_timer.start();
        solver_.step(timer.currentStepLength(), state);
        solver_timer.stop();
        const double st = solver_timer.secsSinceStart();
        std::cout << "Lattice Boltzmann Simulator took: " << st << " seconds." << std::endl;
        stime += st;
        ++timer;
    }
    total_timer.stop();
    std::cout << "Total time taken: " << total_timer.secsSinceStart()
              << "\n Solver Time taken: " << stime;
}
