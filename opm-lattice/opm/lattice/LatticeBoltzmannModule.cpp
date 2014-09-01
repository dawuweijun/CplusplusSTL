#include <opm/lattice/LatticeBoltzmannModule.hpp>
#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>

LatticeBoltzmannModule::LatticeBoltzmannModule()
       : name_("D3Q19")
       , number_of_direction_(19)
       , weight_(std::vector<double>(19))
       , xVelocity_(std::vector<double>(19))
       , yVelocity_(std::vector<double>(19))
       , zVelocity_(std::vector<double>(19))
{
    const double w0 = 1.0 / 3.0;
    const double w1 = 1.0 / 36.0;
    const double w2 = 1.0 / 18.0;
    const double c0 = 0.0;
    const double c2 = 1.0;
    /* Default directional velocities */
    const double c_x[19] = {c2,-c2, c0, c0, c0, c0, c2,-c2, c2,-c2, c2,-c2, c2,-c2, c0, c0, c0, c0, c0};
    const double c_y[19] = {c0, c0, c2,-c2, c0, c0, c2,-c2,-c2, c2, c0, c0, c0, c0, c2,-c2, c2,-c2, c0}; 
    const double c_z[19] = {c0, c0, c0, c0, c2,-c2, c0, c0, c0, c0, c2,-c2,-c2, c2, c2,-c2,-c2, c2, c0}; 
    /* Weight array */
    const double w[19] = {w2, w2, w2, w2, w2, w2, w1, w1, w1, w1, w1, w1, w1, w1, w1, w1, w1, w1, w0};
    xVelocity_ = std::vector<double>(c_x, c_x + 19);
    yVelocity_ = std::vector<double>(c_y, c_y + 19);
    zVelocity_ = std::vector<double>(c_z, c_z + 19);
    weight_    = std::vector<double>(w  , w   + 19);
}
