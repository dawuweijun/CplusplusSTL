#include <utlity>
#include <iostream>

//using std::rel_ops defines the other comparison
//automatically.
class X {
    public:
        bool operator== (const X& x) const;
        bool operator< (const X& x) const;
};

void foo()
{
    using namespace std::rel_ops;
    X x1, x2;
    if (x1 != x2) {
        std:cout << "not equal Works!" << std::endl;
    }
    if (x1 > x2) {
        std::cout << "greater works!" << std::endl;
    }
}
