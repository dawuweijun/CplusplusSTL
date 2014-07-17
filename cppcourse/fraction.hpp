#ifndef FRACTION_INCLUDED
#define FRACTION_INCLUDED
class fraction {
public:
    fraction(int num, int den);
    int num() const;
    int den() const;
    ~fraction();
    fraction add(const fraction& x, const fraction& y);
    fraction reduced();
private:
    int gcd(int m, int n) const;
    int num_;
    int den_;
};
#endif
