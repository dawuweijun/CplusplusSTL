#include <set>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <string>

template<typename T>
inline void print(const T& coll, const std::string& optstr="")
{
    std::cout << optstr;
    for(const auto& elem : coll) {
        std::cout << elem << ' ';
    }
    std::cout << std::endl;
}

int square(int value)
{
    return value * value;
}

int main()
{
    std::set<int> coll1;
    std::vector<int> coll2;

    for (int i = 1; i <= 9; ++i) {
        coll1.insert(i);
    }
    print(coll1, "initialized: ");
    std::transform(coll1.cbegin(), coll1.cend(), std::back_inserter(coll2),square);
    print(coll2, "squared: ");
}
