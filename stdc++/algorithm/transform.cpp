#include <iostream>		//std::cout
#include <algorithm>	//std::tranhsform
#include <vector>		//std::vector
#include <functional>	//std::plus

int op_increase (int i) {return ++i;}

int main()
{
	std::vector<int> foo;
	std::vector<int> bar;
	
	//set some values:
	for (int i=1; i<6; i++){
		foo.push_back(i*10);
	}

	bar.resize(foo.size());

	std::transform(foo.begin(), foo.end(), bar.begin(), op_increase);//bar: 11 21 31 41 51

	std::transform(foo.begin(), foo.end(), bar.begin(), foo.begin(), std::plus<int>());
	
	std::cout << "bar contains: ";
	for (std::vector<int>::iterator it=bar.begin(); it!= bar.end(); ++it){
		std::cout << ' ' << *it;
	}
	std::cout << '\n';

	
	std::cout << "foo contains: ";
	for (std::vector<int>::iterator it=foo.begin(); it!= foo.end(); ++it){
		std::cout << ' ' << *it;
	}
	std::cout << '\n';

	return 0;

}
