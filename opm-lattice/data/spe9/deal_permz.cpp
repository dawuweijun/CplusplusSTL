#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

void read_write(std::ifstream& in, std::ofstream& out)
{
	std::vector<double> perm;
	
	std::istringstream istr;
	string line;
	double tmp;
	out << "PERMZ" << std::endl;
	while (getline(in, line)) {
		istr.str(line);
		while (istr >> tmp){
			perm.push_back(tmp);
		}
		istr.clear();
		line.clear();
	}
	in.close();
	int j = 1, k = 1;
		out << "--LAYER " << k << '\n';
		out << "--ROW   " << j << '\n';
	for (std::vector<int>::size_type i = 0; i != perm.size(); ++i){
		out << std::setprecision(10) << '\t' << perm[i] * 0.01 <<'\t';
		if ((i+1) % 24 == 0){
			out << '\n';
			if ((i+1) % (24 * 25) == 0){
				k++;
				if (k <= 15) 
					out << "--LAYER " << k << '\n';
				else{
					k = 1;
					out << "--LAYER " << k << '\n';
				}
			}
			j++;
			if (j <= 25)
				out << "--ROW  " << j << '\n';
			else{
				j = 1;
				out << "--ROW  " << j << '\n';
			}
		}
		if (((i+1) % 24) % 5 == 0)
			out << '\n';
	}
	out.close();
}
int main (void)
{
	std::ifstream in;
	std::ofstream out;

	in.open ("tmp.dat");
	out.open ("permz.inc");
	if (!in && !out) {
		std::cerr << "error: unable to open input and output file: "
				  << in << out << std::endl;
		return -1;
	}
	read_write(in, out);	
	
/*	if (!infile) {
		return -1;
	}
	while (infile.good() && !infile.eof()) {		
//		getline(infile, line);
		infile.read(perm, sizeof(perm));
		cout << perm;
 		if (line[0] != '-') {
			if (line.find("PERM") != 0){
				for (string::iterator iter = line.begin(); iter != line.end(); ++iter){
					*iter *= 0.1; 
						outfile << line << std::endl;
					}
			}
		}
	}

 	infile.close();
	outfile.close();
*/
 	return 0;
	
}
