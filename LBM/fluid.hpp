#ifndef FLUID_HEADER_INCLUDED
#define FLUID_HEADER_INCLUDED
/*
author: Ming Liu
email:  qilicun@outlook.com
*/
#include "grid.hpp"
#include <vector>
class Grid;
class Fluid{
public:
    Fluid();
    Fluid(Grid& gird, double, double, double, double);
    ~Fluid();
    void init();
    double rho() const;
    double tau() const;
    double mu() const;
    double velmax() const;
    std::vecotr<double> distribution() const;
private:
    double rho_;
    double tau_;
    double mu_;
    double velmax_;
    std::vector<double> velocity_;
    std::vector<std::vector<double> >distribution_;
    Grid& grid_;
    Module& module_;
     
    double feq(double rho, double )
};
#endif //FLUID_HEADER_INCLUDED
