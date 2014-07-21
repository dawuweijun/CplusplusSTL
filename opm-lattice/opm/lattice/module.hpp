#ifndef MODULE_HEADER_INCLUDED
#define MODULE_HEADER_INCLUDED

#include <vector>
#include <string>
#include <cmath>

class Module{
public:
    Module();
    ~Module();
#if 0
    void init()
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
        const double w[19] = {w2, w2, w2, w2, w2, w2, w1, w1, w1, w1, w1, w1, w1, w1, w1, w1, w1, w1, w0};
        for (int i = 0; i < xVelocity_.size(); ++i) {
            xVelocity_[i] = c_x[i];
            yVelocity_[i] = c_y[i];
            zVelocity_[i] = c_z[i];
            weight_[i] = w[i];
        }
    }
#endif    
    inline int numDirection() const 
    {   
        return number_of_direction_;
    }
    std::string name() const {return name_;}
    std::vector<double>& xVelocity() {return xVelocity_;}
    std::vector<double>& yVelocity() {return yVelocity_;}
    std::vector<double>& zVelocity() {return zVelocity_;}
    std::vector<double>& weight() {return weight_;}
   
    const std::vector<double>& xVelocity() const  {return xVelocity_;}
    const std::vector<double>& yVelocity() const {return yVelocity_;}
    const std::vector<double>& zVelocity() const {return zVelocity_;}
    const std::vector<double>& weight() const {return weight_;}
private:
    std::string name_;
    int number_of_direction_; 
    std::vector<double> weight_;
    std::vector<double> xVelocity_;
    std::vector<double> yVelocity_;
    std::vector<double> zVelocity_;
};
#endif//MODULE_HEADER_INCLUDED
