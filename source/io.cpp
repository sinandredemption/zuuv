#include "io.h"
#include <cstdio>
#include <iostream>
#include <fstream>

void read_from_file(mpq_class* frac, std::string file_name, bool dec)
{
	std::cout << "Reading from file: " << file_name << "...";
	std::ifstream file(file_name);

	if (!file)
		throw "read_from_file(): FAILED (couldn't open file)";
	
	std::string str;
	getline(file, str);

	if (str[1] == '.') {
		str[1] = str[0];
		str[0] = '0';
	}

	mpz_set_str(frac->get_num_mpz_t(), str.c_str(), dec ? 10 : 16);

	file.close();

	std::cout << " dec digits: \t" << mpz_sizeinbase(frac->get_num_mpz_t(), 10)
		<< " | hex digits: \t" << mpz_sizeinbase(frac->get_num_mpz_t(), 16) << std::endl;

	frac->get_den() = 1;
	frac->get_den() <<= mpz_sizeinbase(frac->get_num_mpz_t(), 2) - 2;

	std::cout << "Reading completed. Fraction ~= " << frac->get_d() << std::endl;
}

void output_cf_terms_list(const CFTerms& cf_term_list, std::string file_name, bool verbose) {
	if (verbose)
		std::cout << "Writing " << cf_term_list.size() << " terms to file: " << file_name << "...";

	std::ofstream out(std::string("iterations/") + file_name);

	if (!out) {
		if (verbose)
		std::cerr << "FAILED (couldn't open file)\n";
		throw "Couldn't open file";
		return;
	}

	for (auto c : cf_term_list)
		out << c << '\n';

	if (verbose)
		std::cout << "COMPLETED" << std::endl;
	out.close();
}
