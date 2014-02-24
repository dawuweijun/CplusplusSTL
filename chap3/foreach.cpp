#include <iostream>
#include <vector>

int main()
{
	for (int i : {2, 3, 4, 5, 6, 7}) {
		std::cout << i << std::endl;
	}

	std::vector<double> vec {2,3,4,5,6,7,8};
	for (auto& elem : vec) {
		elem *= 3;
	}
	
	for (auto& elem : vec) {
		std::cout << elem << std::endl;
	}
	
	int array[] = {1, 2, 3, 4, 5};
	long sum = 0;
	for (int x : array) {
		sum += x;
	}

	for (auto elem : {sum, sum*2, sum*4}) {
		std::cout << elem << std::endl;
	}
}
