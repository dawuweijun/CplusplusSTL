#include <iostream>

class P
{
	public:
		P(int, int);
		P(std::initializer_list<int>);
};

class Q
{
	public:
		Q(int a, int b);
		explicit Q(int a, int b, int c);
}

int main()
{
	P p(77,5);
	P q{77,5};
	P r{77,5,42};
	P p(77,5);
	P s = {77, 5};

	Q x(77, 5);
	Q y{77, 5};
	Q z {77, 5, 42};
	Q v = {77, 5};
	Q w = {77, 5, 42};
}
