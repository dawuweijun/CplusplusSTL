#include <iostream>
#include <ctime>
using namespace std;
int main()
{
    clock_t start, finish;
    start = clock();
    std::cout << "hahhaah: "<< std::endl;
    finish = clock();
    std::cout << "time: " << finish-start << std::endl;
}
