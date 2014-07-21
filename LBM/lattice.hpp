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
    Lattice(const Grid& grid, const Module& module, const Fluid& red, const Fluid& blue);
    ~Lattice();
    void 
    step(const Grid& grid, Fluid& red, Fluid& blue);
private:
    double externalForce_;
    double fluxForce_;
    double flux_;
    std::vector<double> gff_;
    std::vector<double> gfs_;
    Grid& grid_;
    Module& module_;
    Fluid& red_;
    Fluid& blue_;
    void potential();
    void propagationBySwap(const Grid& grid, Fluid& fluid);
    void collision(const Grid& grid, Fluid& red, Fluid& blue);
    void streamingSwap(const Grid& grid, Fluid& red, Fluid& blue);
    void massMomentumCalc(const std::vector<Fluid>& fluid, const Grid& grid);
    void pressureCalc(const std::vector<Fluid>& fluid, const Grid& grid);
}
#endif //LATTICE_HEADER_INCLUDED
