#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

int read_write(std::istream& is, std::ostream& out)
{
	std::vector<double> zcorn;
	zcorn.clear();
	while (is) {
		double line;
		is >> line;
		if (is.eof()) 
			break;
		zcorn.push_back(line);
	}
	std::cout << zcorn.size() << std::endl;
	std::vector<int>::size_type ix;
	for (ix = 0; ix != zcorn.size(); ++ix) {	
		out << zcorn[ix] << ' ';
		if ((ix + 1) % 6 == 0)
			out << std::endl; 
	}
	out << std::endl;
}

int main(void)
{
	std::ifstream in;
	std::ofstream out;
	std::string read_file="test";
	std::string write_file="tmp";
	in.open(read_file.c_str());
	out.open(write_file.c_str());
	if (!in && !out) {
		std::cerr << "error: unable to open input and output file: "
				  << in << out << std::endl;
		return -1;
	}
	read_write(in, out);

	return 0;	
}
