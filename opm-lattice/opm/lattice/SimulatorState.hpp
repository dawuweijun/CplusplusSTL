#ifndef SIMULATOR_STATE_HEADER_INCLUDED
#define SIMULATOR_STATE_HEADER_INCLUDED
#include <vector>
#include <opm/lattice/GridManager.hpp>
class SimulatorState
{
public:
    void init (const GridManager& grid, const LatticeBoltzmannModule& module);
    int numPhases() const { return num_phases_; }
    std::vector<double>& pressure() { return press_; }
    std::vector<double>& velocity() { return vel_; }
    std::vector<std::vector<double>>& distribution { return distr_; }
    const std::vector<double>& pressure() const { return press_; }
    const std::vector<double>& velocity() const { return vel_; }
    const std::vector<std::vector<double>>& distribution const { return distr_; }
private:
    int num_phases_;
    std::vector<double> press_;
    std::vector<double> vel_;
    std::vector<std::vector<double>> distr_;
};

#endif // SIMULATOR_STATE_HEADER_INCLUDED
