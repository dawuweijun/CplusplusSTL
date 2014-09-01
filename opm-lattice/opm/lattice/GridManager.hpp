#ifndef GRIDMANAGER_HEADER_INCLUDED
#define GRIDMANAGER_HEADER_INCLUDED
/*
author: Ming Liu
email:  qilicun@outlook.com
*/
#include <vector>

class GridManager{
public:
    GridManager(const int NX, const int NY, const int NZ);
    GridManager(const int NX, const int NY);
    const int dimension() const;
    const int spaceDim() const;
    const int index(const int x, const int y, const int z) const;

    const std::vector<int>& boundary() const { return boundary_; }
    std::vector<int>& boundary() {return boundary_; }
    const int NX() const { return NX_; }
    const int NY() const { return NY_; }
    const int NZ() const { return NZ_; }
     
private:
    const int NX_;
    const int NY_;
    const int NZ_;
    std::vector<int> boundary_;
};
#endif //GRID_HEADER_INCLUDED
