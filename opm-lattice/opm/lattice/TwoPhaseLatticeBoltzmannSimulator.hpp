#ifndef TWOPHASELATTICEBOLTZMANNSIMULATOR_HEADER_INCLUDED
#define TWOPHASELATTICEBOLTZMANNSIMULATOR_HEADER_INCLUDED

#include <opm/lattice/FluidProperties.hpp>
#include <opm/lattice/GridManager.hpp>
#include <opm/lattice/LatticeBoltzmannModule.hpp>
#include <opm/lattice/LatticeBoltzmannSolver.hpp>
#include <opm/lattice/SimulatorState.hpp>
#include <opm/lattice/writeVtkData.hpp>
#include <opm/lattice/DataMap.hpp>
#include <opm/lattice/SimulatorTimer.hpp>
#include <boost/filesystem.hpp>
#include <opm/lattice/StopWatch.hpp>
#include <iostream>
#include <string>
class TwoPhaseLatticeBoltzmannSimulator 
{
public:
    TwoPhaseLatticeBoltzmannSimulator(const GridManager& grid,
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
    LatticeBoltzmannSolver solver_;

    void outputStateVtk(const GridManager& grid,
                        const SimulatorState& state,
                        const int step,
                        const std::string& output_dir);
    void outputStateMatlab(const GridManager& grid,
                           const SimulatorState& state,
                           const int step,
                           const std::string& output_dir);
};
#endif //TWOPHASELATTICEBOLTZMANNSIMULATOR_HEADER_INCLUDED
