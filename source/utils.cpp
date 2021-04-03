#include "utils.h"
#include <sstream>
#include <fstream>
#include <iterator>

double wall_clock() {
	//  Get the clock in seconds.
	auto ratio_object = std::chrono::high_resolution_clock::period();
	double ratio = (double)ratio_object.num / ratio_object.den;
	return std::chrono::high_resolution_clock::now().time_since_epoch().count() * ratio * 1000.;
}

mpz_class mpz_rand(size_t bytes) {
	mpz_class out;
	size_t length = bytes / 8;
	auto ptr = mpz_limbs_write(out.get_mpz_t(), length);
	for (int i = 0; i < length; ++i) {
		unsigned long long rand64 = rand();
		rand64 = (rand64 << 15) | rand();
		rand64 = (rand64 << 15) | rand();
		rand64 = (rand64 << 15) | rand();
		rand64 = (rand64 << 15) | rand();
		ptr[i] = rand64;
	}
	mpz_limbs_finish(out.get_mpz_t(), length);
	return out;
}

unsigned count_lines(std::string filename)
{
	std::ifstream myfile(filename);

	// new lines will be skipped unless we stop it from happening:    
	myfile.unsetf(std::ios_base::skipws);

	// count the newlines with an algorithm specialized for counting:
	unsigned line_count = std::count(
		std::istream_iterator<char>(myfile),
		std::istream_iterator<char>(),
		'\n');

	return line_count;
}

bool file_exists(std::string filename)
{
	std::ifstream file(filename);
	return file.good();
}

