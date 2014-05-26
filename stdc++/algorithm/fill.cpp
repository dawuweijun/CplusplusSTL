#include <iostream> 	//std::cout
#include <vector>		//std::vector
#include <algorithm>	//std::fill

int main()
{
	std::vector<int> myvector (8);

	std::fill (myvector.begin(), myvector.begin()+4, 5);
	std::fill (myvector.begin()+3, myvector.end()-2, 8);
	
	std::cout << "myvector contains: ";
	for (std::vector<int>::iterator it=myvector.begin(); it!=myvector.end(); ++it){
		std::cout << ' ' << *it;
	}
	std::cout << '\n';
	
	return 0;
}

