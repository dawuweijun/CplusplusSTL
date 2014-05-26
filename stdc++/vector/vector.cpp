#include <iostream>
#include <vector>

void print (const std::vector<double>& a)
{
	std::cout << "Now the vector is:" << std::endl;
	std::vector<double>::const_iterator iter;
	for (iter = a.begin(); iter != a.end(); ++iter)
		std::cout << *iter << "  ";
	std::cout << std::endl;
}
int main()
{
	std::vector<double> a;
	for (int i = 0; i < 4; i++)
		a.push_back(i);
	std::cout << "init vector" << std::endl;
	print(a);	
	for (int j = 0; j < 4; j++)
		a.push_back(j+4);
	std::cout << "after re-push_back" << std::endl;
	print(a);
} 	
