#ifndef MODULE_HEADER_INCLUDED
#define MODULE_HEADER_INCLUDED
#include <opm/lattice/utility/SimulatorState.hpp>

#include <vector>
#include <string>

class LatticeBoltzmannModule{
public:
    LatticeBoltzmannModule();
    const int numDirection() const { return number_of_direction_; }
    const std::string name() const { return name_; }
    std::vector<double>& xVelocity() { return xVelocity_; }
    std::vector<double>& yVelocity() { return yVelocity_; }
    std::vector<double>& zVelocity() { return zVelocity_; }
    std::vector<double>& weight() { return weight_; }
   
    const std::vector<double>& xVelocity() const  { return xVelocity_; }
    const std::vector<double>& yVelocity() const { return yVelocity_; }
    const std::vector<double>& zVelocity() const { return zVelocity_; }
    const std::vector<double>& weight() const { return weight_; }

private:
    std::string name_;
    int number_of_direction_; 
    std::vector<double> weight_;
    std::vector<double> xVelocity_;
    std::vector<double> yVelocity_;
    std::vector<double> zVelocity_;
};
#endif//MODULE_HEADER_INCLUDED
