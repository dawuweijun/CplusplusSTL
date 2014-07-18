#ifndef FLUID_HEADER_INCLUDED
#define FLUID_HEADER_INCLUDED
/*
author: Ming Liu
email:  qilicun@outlook.com
*/
#include "grid.hpp"
#include "module.hpp"
#include <vector>
class Grid;
class Module;
class Fluid{
public:
    Fluid(const Grid& gird, const Module& module, const double, const double, const double, const double);
    ~Fluid();
    enum Location = {Inside = 0, Outside = 1};
    void init(const double x1, const double x2, const Location loc);
    double rho() const;
    double tau() const;
    double mu() const;
    double velmax() const;
    inline double lambda() const {return -2.0 / (6*tau_*mu_ + 1);}
    const std::vecotr<std::vector<double> >& distribution() const {return distribution_;}
    std::vecotr<std::vector<double> >& distribution() { return distribution_;}
private:
    double rho_;
    double tau_;
    double mu_;
    double velmax_;
    std::vector<double> velocity_;
    std::vector<double> density_;
    std::vector<std::vector<double> >distribution_;
    Grid& grid_;
    Module& module_;
    std::vector<double> feq();
};
#endif //FLUID_HEADER_INCLUDED
