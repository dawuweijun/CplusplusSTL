#include <iostream>

const int max(const int a, const int b)
{
    if (a >= b) {
        return a;
    } else {
        return b;
    }
}
int main()
{
    int a, b ,c;
    std::cout << "Please input 3 numbers: ";
    std::cin >> a;
    std::cin >> b;
    std::cin >> c;
    int max_num;
    max_num = max(a, b);
    max_num = max(max_num, c);
    std::cout << "The maxmium of " << a <<","<< b << "," << c << " is: " << max_num << std::endl;
    return 0;
}
