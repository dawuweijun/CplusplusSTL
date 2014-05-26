#include <iostream>
#include <cstring>
#include <string>

int main()
{
	std::string str ("Please split this sectence into tokens");
	
	char *cstr = new char [str.length() + 1];
	std::strcpy (cstr, str.c_str());
	
	char *p = std::strtok (cstr, " ");
	while (p != NULL){
		std::cout << p << '\n';
		p= strtok(NULL, " ");
	}

	delete[] cstr;
	return 0;

}
