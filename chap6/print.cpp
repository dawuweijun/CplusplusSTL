#include <iostream>
#include <string>
#include <vector>
template<typename T>
inline void print(const T& coll, const std::string& optstr="")
{
    std::cout << optstr;
    for(const auto& elem : coll) {
        std::cout << elem << ' ';
    }
    std::cout << std::endl;
}

int main()
{
    std::vector<int> coll{1,2,3,4,5};
    
    print(coll, "All the ints are:\n");
}
