#include <iostream>
#include <vector>

std::vector<int> fibonacci(int n)
{
    std::vector<int> fibo;
    int a = 0;
    int b = 1;
    while (a < n) {
        std::cout << a << " " << b << std::endl;
        fibo.push_back(a);
        b = a + b;
        a = b;
    }
    return fibo;
}

int fibonacci2(int n)
{
    if (n == 0) {
        std::cout << "number must be greater 0.\n";
        exit(1);
    }
    if (n == 1 || n == 2)
        return 1;
    return fibonacci2(n-1) + fibonacci2(n-2);
}
int main()
{
    std::cout << "Please input a number: ";
    int n;
    std::cin >> n;
    std::vector<int> fibo;
    
    fibo = fibonacci(n);
    for (auto i : fibo) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
    int f = fibonacci2(n);
    std::cout << f << " ";
    std::cout << std::endl;

    return 0;
}
