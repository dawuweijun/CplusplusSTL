#ifndef GRID_HEADER_INCLUDED
#define GRID_HEADER_INCLUDED
/*
author: Ming Liu
email:  qilicun@outlook.com
*/
#include <vector>
class Grid{
public:
    Grid(double NX, double NY, double NZ);
    Grid(double NX, double NY);
    ~Grid();
    double dimension() const;
    double spaceDim() const;
    double index(double x, double y, double z) const;
    std::vector<int> boundary();
private:
    double NX_;
    double NY_;
    double NZ_;
    std::vector<int> boundary_;
};
#endif //GRID_HEADER_INCLUDED
