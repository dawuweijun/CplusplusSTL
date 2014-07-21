#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

int main(void)
{
	const char *filename = "poro.inc";	
	double p[15]={0.087, 0.097, 0.111, 0.16, 0.13, 0.17, 0.17, 0.08, 0.14, 0.13,0.12,0.105,0.12,0.116,0.157};	
	
	std::ofstream out;
	out.open(filename);
	if (!out) {
		std::cerr << "error: unable to open input and output file: "
				  << out << std::endl;
		return -1;
	}

	out << "PORO" << std::endl;
	for (int k = 0; k < 15; k++){
		out << "--LAYER " << k + 1<< '\n';
		for (int j = 0; j < 25; j++){
			out << "--ROW   " << j + 1<< '\n';
			for (int i = 0; i < 24; i++){
				out << p[k] << '\t';
				if ((i+1) % 5 == 0)
				out << '\n';	 
			}
			out << '\n';
		}
	}
	out.close();
 	return 0;
}
