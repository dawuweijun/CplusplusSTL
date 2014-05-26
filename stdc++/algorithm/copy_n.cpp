#include <iostream>
#include <algorithm>
#include <vector>

int main()
{
    int myints[] = {10,20,30,50,60};
    std::vector<int> myvector;
    
//    myvector.resize(5);
    
    std::copy_n(myints, 9, myvector.begin());
    std::cout << "myvector:\n";
    
for (std::vector<int>::iterator it = myvector.begin(); it!=myvector.end(); ++it)
    std::cout << ' ' << *it;

    std::cout << '\n';
    return 0;
    
}
