#include <vector>
#include <iostream>

int main(void)
{
	int num = 5;
	std::vector<int> a;
	for (int i = 1; i <= num; ++i)
		a.push_back(i);
	std::vector<int>::const_iterator it;

	std::cout << "Before resize the vector:\n";

	for (it = a.begin(); it != a.end(); ++it) {
		std::cout<< *it << " ";
	}
	std::cout << std::endl;
	std::cout << "After resize the vector:\n";
	int re = 8;
	a.clear();
	a.resize(re, 100);
	for (it = a.begin(); it != a.end(); ++it) {
		std::cout<< *it << " ";
	}
	std::cout << std::endl;	
}
