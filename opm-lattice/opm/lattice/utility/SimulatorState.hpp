#ifndef SIMULATOR_STATE_HEADER_INCLUDED
#define SIMULATOR_STATE_HEADER_INCLUDED

#include <vector>
class GridManager;
class LatticeBoltzmannModule;

class SimulatorState
{
public:
    void init (const GridManager& grid, const LatticeBoltzmannModule& module);
    int numPhases() const { return num_phases_; }
    std::vector<double>& pressure() { return press_; }
    std::vector<double>& velocity() { return vel_; }
    std::vector<double>& redDistr() { return red_distr_;}
    std::vector<double>& blueDistr() { return blue_distr_;}
    const std::vector<double>& pressure() const { return press_; }
    const std::vector<double>& velocity() const { return vel_; }
    const std::vector<double>& redDistr() const { return red_distr_;}
    const std::vector<double>& blueDistr() const { return blue_distr_;}
private:
    int num_phases_;
    std::vector<double> press_;
    std::vector<double> vel_;
    std::vector<double> red_distr_;
    std::vector<double> blue_distr_;
};

#endif // SIMULATOR_STATE_HEADER_INCLUDED
