#include <vector>
#include <iostream>

using namespace std;

int main()
{
	std::vector<int> first;
	std::vector<int> second;
	std::vector<int> third;
	
	first.assign(7,100);
	
	std::vector<int>::iterator it;
	it = first.begin()+1;

	second.assign(it, first.end()-1);

	int myints[]= {1774, 7, 4};
	third.assign (myints, myints+3);

	std::cout << "Size of first: " << int (first.size()) << '\n';
	std::cout << "Value of first: " << '\n';
	for (it = first.begin(); it!=first.end(); ++it){
		std::cout << *it << " ";
	}
	std::cout << '\n';
	std::cout << "Size of second: " << int (second.size()) << '\n';
	std::cout << "Value of second: " << '\n';
	for (it = second.begin(); it!=second.end(); ++it){
		std::cout << *it << " ";
	}

	std::cout << '\n';
	std::cout << "Size of third: " << int (third.size()) << '\n';
	std::cout << "Value of third: " << '\n';
	for (it = third.begin(); it!=third.end(); ++it){
		std::cout << *it << " ";
	}
	
	return 0;
}
