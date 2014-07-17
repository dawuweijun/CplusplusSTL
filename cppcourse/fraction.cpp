#include "fraction.hpp"
#include <cmath>
#include <iostream>

fraction::fraction(int num, int den)
: num_(num)
, den_(den)
{
    if (den == 0) {
        std::cout <<"invalide den:\n";
        exit(1);
    }
}

fraction::~fraction()
{
}

int
fraction::num() const
{
    return num_;
}

int 
fraction::den() const
{
    return den_;
}

fraction fraction::add(const fraction& x, const fraction& y)
{
    num_ = x.num() * y.den() + x.den() * y.num();
    den_ = x.den() * y.den();
    return reduced();
}

fraction fraction::reduced()
{
   int num = gcd(num_, den_);
   num_ /= num;
   den_ /= num;
   return fraction(num_, den_);
}

int
fraction::gcd(int m, int n) const
{
    int max = std::max(m, n);
    int min = std::min(m, n);
    if (max % min == 0) {
        return min;
    } else {
        return gcd(min, max % min);
    }
}

int main()
{
    int num, den;
    std::cout << "Please intput two fractions: ";
    std::cout << "Frist(num den): ";
    std::cin >> num;
    std::cin >> den;
    fraction x(num, den);
    std::cout << "Second(num den): ";
    std::cin >> num;
    std::cin >> den;
    fraction y(num, den);
    fraction z(1,1);
    z.add(x, y);
    std::cout << "Sum of them is: " << z.num() << "/" << z.den() << std::endl;

    return 0;
}
