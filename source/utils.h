#ifndef INC_UTILS_H
#define INC_UTILS_H
#include <chrono>
#include <mpirxx.h>
#include <sstream>
#include <string>

double wall_clock();
mpz_class mpz_rand(size_t length);
template <typename T>
inline std::string to_str(T t) {
	std::stringstream os;
	os << t;
	return os.str();
}
unsigned count_lines(std::string filename);
bool file_exists(std::string filename);

#endif
