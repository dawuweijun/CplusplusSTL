#ifndef MODULE_HEADER_INCLUDED
#define MODULE_HEADER_INCLUDED

#include <vector>
#include <string>
class Module{
public:
    Module();
    ~Module();
    int numDirection() const;
  //  std::string& name() const;
    std::vector<double> xVelocity() const;
    std::vector<double> yVelocity() const;
    std::vector<double> zVelocity() const;
    std::vector<double> weight() const;
   
private:
//    std::string name_;
    int number_of_direction_; 
    std::vector<double> weight_;
    std::vector<double> xVelocity_;
    std::vector<double> yVelocity_;
    std::vector<double> zVelocity_;
};
#endif//MODULE_HEADER_INCLUDED
