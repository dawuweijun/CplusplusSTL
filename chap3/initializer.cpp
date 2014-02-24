#include <vector>
#include <iostream>

void print (std::initializer_list<int> vals)
{
	for (auto p=vals.begin(); p!=vals.end(); ++p) {
		std::cout << *p << "\n";
	}
}

int main()
{
	print({12,3,5,7,11,13,17});
	return 0;
}
