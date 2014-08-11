#ifndef TWOPHASELATTICEBOLTZMANNSIMULATOR_HEADER_INCLUDED
#define TWOPHASELATTICEBOLTZMANNSIMULATOR_HEADER_INCLUDED

#include <memory>
#include <string>

class FluidProperties;
class GridManager;
class LatticeBoltzmannModule;
class LatticeBoltzmannSolver;
class SimulatorState;
class SimulatorTimer;

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
    class Impl;
    std::shared_ptr<Impl> pimpl_;

};
#endif //TWOPHASELATTICEBOLTZMANNSIMULATOR_HEADER_INCLUDED
