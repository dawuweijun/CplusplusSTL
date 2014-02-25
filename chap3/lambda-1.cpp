#include <iostream>

int main()
{
	int x = 0;
	int y = 42;
	auto qqq = [x, &y] {
					std::cout << "x: " << x << std::endl;
					std::cout << "y: " << y << std::endl;
	};
	x = y = 77;
	qqq();
	qqq();
	std::cout << "final y: " << y << std::endl;


	int id = 0;
	auto f = [id] () mutable {
						std::cout << "id: " << std::endl;
						++id;
	};
	id = 42;
	f();
	f();
	f();
	std::cout << "final id: " << id << std::endl;
	return 0;
}
