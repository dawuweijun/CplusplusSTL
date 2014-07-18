#ifndef LATTICE_HEADER_INCLUDED
#define LATTICE_HEADER_INCLUDED

#include "grid.hpp"
#include "module.hpp"
#include "fluid.hpp"
#include <vector>

class Grid;
class Module;
class Fluid;
class Lattice {
public:
    Lattice(const Grid& grid, const Module& module, const std::vector<Fluid>& fluid);
    ~Lattice();
//    void massMomentumCalc(const std::vector<Fluid>& fluid, const Grid& grid);
//    void pressureCalc(const std::vector<Fluid>& fluid, const Grid& grid);
    void collision(const Grid& grid, const std::vector<Fluid>& fluid);
    void streamingSwap(const Grid& grid, const std::vector<Fluid>& fluid);
private:
    double externalForce_;
    double fluxForce_;
    double flux_;
    std::vector<double> gff_;
    std::vector<double> gfs_;
    Grid& grid_;
    Module& module_;
    std::vector<Fluid>& fluid_;
    void potential();
    void propagationBySwap(const Grid& grid, Fluid& fluid);
}
#endif //LATTICE_HEADER_INCLUDED
