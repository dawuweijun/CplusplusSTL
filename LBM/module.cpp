#include "module.hpp"
#include <string>
#include <cmath>
#include <vector>
#include <iostream>

Module::Module()
//       :name_("D3Q19")
       :number_of_direction_(19)
       ,weight_(19)
       ,xVelocity_(19)
       ,yVelocity_(19)
       ,zVelocity_(19)
{
   const double w0 = 1.0 / 3.0;
   const double w1 = 1.0 / 36.0;
   const double w2 = 1.0 / 18.0;
   const double c0 = 0.0;
   const double c1 = std::sqrt(2);
   const double c2 = 1.0;
    /* Default directional velocities */
   const double c_x[19] = {c2,-c2, c0, c0, c0, c0, c2,-c2, c2,-c2, c2,-c2, c2,-c2, c0, c0, c0, c0, c0};
   const double c_y[19] = {c0, c0, c2,-c2, c0, c0, c2,-c2,-c2, c2, c0, c0, c0, c0, c2,-c2, c2,-c2, c0}; 
   const double c_z[19] = {c0, c0, c0, c0, c2,-c2, c0, c0, c0, c0, c2,-c2,-c2, c2, c2,-c2,-c2, c2, c0}; 
    /* Weight array */
    const double w[19] = {w2, w2, w2, w2, w2, w2,                          // Sublattice 1
			    w1, w1, w1, w1, w1, w1, w1, w1, w1, w1, w1, w1,  // Sublattice 2
			    w0};                                             // Sublattice 0 (rest particle)
     for (int i = 0; i < xVelocity_.size(); ++i) {
        xVelocity_[i] = c_x[i];
        yVelocity_[i] = c_y[i];
        zVelocity_[i] = c_z[i];
        weight_[i] = w[i];
     }
}
Module::~Module()
{}

int
Module::numDirection() const
{
    return number_of_direction_;
}
/*
std::string& Module::name() const
{
    return name_;
}
*/
std::vector<double> Module::xVelocity() const 
{
    return xVelocity_;
}
std::vector<double> Module::yVelocity() const 
{
    return yVelocity_;
}
std::vector<double> Module::zVelocity() const 
{
    return zVelocity_;
}
std::vector<double> Module::weight() const 
{
    return weight_;
}

void print (const std::string s, std::vector<double>& vec) 
{
    std::cout << s << std::endl;
    for (auto elem : vec) {
        std::cout << elem << " " << std::endl;
    }
}
int main()
{
    Module module;
    const int num = module.numDirection();
    std::cout << "direction size: " << num<< std::endl;
//    std::cout << module.name() << std::endl;
    std::vector<double> xvel = module.xVelocity();
    std::vector<double> yvel = module.yVelocity();
    std::vector<double> zvel = module.zVelocity();
    std::vector<double> w = module.weight();
    print("xvel", xvel); 
    print("yvel", yvel); 
    print("zvel", zvel); 
    print("weight",w); 
}
