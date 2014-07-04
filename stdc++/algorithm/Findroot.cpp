#include <iostream>
#include <cmath>

double fun(double x)
{
	double A1 = 4.14898454e05;
	double A2 = 1.25577380e01;
	double x0 = 1;
	double p  = 9.99841781e-01;

	return  A1 + (A2 - A1) / (1 + std::pow(10, (log10(x0)-x)*p)) - 50.;
}

int main()
{
	double x1 = 4;
	double x2 = 5;
	if (fun(x1) * fun(x2) > 1e-8) {
		std::cout << "no root" << std::endl;
	} else {
		while ((x2 - x1) > 1e-8)
		{	
			double mid = (x1 + x2) / 2;
			if (fun(mid) * fun(x2) < 0)
				x1 = mid;
			else
				x2 = mid;
		}
	}
	std::cout << "The root is: " << x1 << " " << x2 << std::endl;
	std::cout << fun(4) + 50 << std::endl;

	return 0;
}
