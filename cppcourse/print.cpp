#include <iostream>
#include <string>

void print(const std::string s, const int n) 
{
    for (int i = 0; i < n; ++i) {
        std::cout << "*";
    }
    std::cout << std::endl;
}
int main ()
{
    const std::string s = "*";
    for (int n = 1; n <10;) {
        print(s, n);
        n = n + 2;
    }
    for (int n = 9; n >0;) {
        print(s, n);
        n = n - 2;
    }
}
