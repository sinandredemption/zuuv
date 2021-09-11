#ifndef INC_UTILS_H
#define INC_UTILS_H
#include <chrono>
#include <mpirxx.h>
#include <sstream>
#include <string>

double wall_clock();
mpz_class mpz_rand(size_t length);
unsigned count_lines(std::string filename);
bool file_exists(std::string filename);
size_t peak_ram_usage(size_t fraction_size);
uint64_t parse_shorthand_num(std::string str);
std::string format_time(double ms);

#endif
