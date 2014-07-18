#include <iostream>

bool isPrime(int n)
{
    int m = std::sqrt(n);
    int i;
    for (i = 2; i < m+1; ++i) {
        if (n % i == 0) {
            break;
        }
    }
    if (i > m) {
        return true;
    } else {
        return false;
    }
}

int main()
{
    std::cout << "Please input a number: ";
    int n;
    std::cin >> n;

    bool flag = isPrime(n);
    if (flag) {
        std::cout << n << " is a prime." << std::endl;
    } else {
        std::cout << n << " is not a prime." << std::endl;
    }

    return 0;
}
