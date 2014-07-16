#include <forward_list>
#include <iostream>
int main()
{
    std::forward_list<int> l1 = {1,2,3,4,5};
    std::forward_list<int> l2 = {97, 98, 99};

    //find 3 in l1
    auto pos1 = l1.before_begin();
    for (auto pb1 = l1.begin(); pb1 != l1.end(); ++pb1, ++pos1) {
        if (*pb1 == 3) {
            break; // found
        }
    }

    //find 99 in l2
    auto pos2 = l2.before_begin();
    for (auto pb2 = l2.begin(); pb2 != l2.end(); ++pb2, ++pos2) {
        if (*pb2 == 99) {
            break; // found
        }
    }
    l1.splice_after(pos2, l2, pos1);
    //print
    for (auto pos = l1.begin(); pos != l1.end(); ++pos) {
        std::cout << *pos << " "; 
    }
    std::cout << std::endl;

    for (auto pos = l2.begin(); pos != l2.end(); ++pos) {
        std::cout << *pos << " ";
    }
    std::cout << std::endl;

    return 0;
}

