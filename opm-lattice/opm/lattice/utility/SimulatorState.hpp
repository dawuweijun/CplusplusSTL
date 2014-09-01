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
    std::vector<double>& redDensity() { return red_density_;}
    std::vector<double>& blueDensity() { return blue_density_;}
    std::vector<double>& redDist() { return red_dist_;}
    std::vector<double>& blueDist() { return blue_dist_;}
    const std::vector<double>& pressure() const { return press_; }
    const std::vector<double>& redDensity() const { return red_density_;}
    const std::vector<double>& blueDensity() const { return blue_density_;}
    const std::vector<double>& redDist() const { return red_dist_;}
    const std::vector<double>& blueDist() const { return blue_dist_;}
private:
    int num_phases_;
    std::vector<double> press_;
    std::vector<double> red_density_;
    std::vector<double> blue_density_;
    std::vector<double> red_dist_;
    std::vector<double> blue_dist_;
};

#endif // SIMULATOR_STATE_HEADER_INCLUDED
