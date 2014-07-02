//
#include <vector>
#include <algorithm>
#include <iostream>

class PrintInt {
public:
    void operator () (int elem) const {
	std::cout << elem << ' ';
    }
};

int main()
{
    std::vector<int> coll;

    for (int i = 1; i <= 9; ++i) {
	coll.push_back(i);
    }

    std::for_each(coll.cbegin(), coll.cend(), PrintInt());
    std::cout << std::endl;
}
