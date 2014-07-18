#ifndef GRID_HEADER_INCLUDED
#define GRID_HEADER_INCLUDED
/*
author: Ming Liu
email:  qilicun@outlook.com
*/
#include <vector>
class Grid{
public:
    Grid(int NX, int NY, int NZ);
    Grid(int NX, int NY);
    ~Grid();
    int dimension() const;
    int spaceDim() const;
    int index(const int x, const int y, const int z) const;
    std::vector<int> boundary();
    const int NX() const { return NX_; }
    const int NY() const { return NY_; }
    const int NZ() const { return NZ_; }
     
private:
    int NX_;
    int NY_;
    int NZ_;
    std::vector<int> boundary_;
};
#endif //GRID_HEADER_INCLUDED
