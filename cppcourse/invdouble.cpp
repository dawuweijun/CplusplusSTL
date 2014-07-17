#include <iostream>
#include <algorithm>
#include <string>
int copy (int i) {return i;}
int inv_double (const int n)
{
    std::string s1, s2;
    s1 = std::to_string(n);
    std::transform(s1.rbegin(), s1.rend(), s2.begin(), copy);
    return 2*std::stoi(s2);
}

int main ()
{
    int n;
    std::cout << "Please input a number: ";
    std::cin >> n;
    int nn = inv_double (n);
    std::cout << "inverse double is: " << nn << std::endl;
    return 0;
}
