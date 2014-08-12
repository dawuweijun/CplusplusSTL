#include <opm/lattice/LatticeBoltzmannModule.hpp>
#include <string>
#include <vector>

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
//   const double c1 = std::sqrt(2);
    const double c2 = 1.0;
    /* Default directional velocities */
    const double c_x[19] = {c2,-c2, c0, c0, c0, c0, c2,-c2, c2,-c2, c2,-c2, c2,-c2, c0, c0, c0, c0, c0};
    const double c_y[19] = {c0, c0, c2,-c2, c0, c0, c2,-c2,-c2, c2, c0, c0, c0, c0, c2,-c2, c2,-c2, c0}; 
    const double c_z[19] = {c0, c0, c0, c0, c2,-c2, c0, c0, c0, c0, c2,-c2,-c2, c2, c2,-c2,-c2, c2, c0}; 
    /* Weight array */
    const double w[19] = {w2, w2, w2, w2, w2, w2,                          // Sublattice 1
               			 w1, w1, w1, w1, w1, w1, w1, w1, w1, w1, w1, w1,  // Sublattice 2
			             w0                                             // Sublattice 0 (rest particle)
                         };
    xVelocity_ = std::vector<double>(c_x, c_x + sizeof(c_x[0]));
    yVelocity_ = std::vector<double>(c_y, c_y + sizeof(c_y[0]));
    zVelocity_ = std::vector<double>(c_z, c_z + sizeof(c_z[0]));
    weight_    = std::vector<double>(w  , w   + sizeof(w[0]));
}


/*
void print (const std::string s, std::vector<double>& vec) 
{
    std::cout << s << std::endl;
    for (auto elem : vec) {
        std::cout << elem << " "; 
    }
    std::cout << std::endl;
}
int main()
{
    LatticeBoltzmannModule module;
    const int num = module.numDirection();
    std::cout << "LatticeBoltzmannModule's name is: " << module.name() << ". It has " << num<< " directions." << std::endl;
    std::vector<double>& xvel = module.xVelocity();
    std::vector<double>& yvel = module.yVelocity();
    std::vector<double>& zvel = module.zVelocity();
    std::vector<double>& w = module.weight();
    print("xvel", xvel); 
    print("yvel", yvel); 
    print("zvel", zvel); 
    print("weight",w); 
}
*/
