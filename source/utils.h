#ifndef INC_UTILS_H
#define INC_UTILS_H
#include <chrono>
#include <mpirxx.h>
#include <sstream>
#include <string>

// High resolution clock (returns time since epoch in seconds)
double wall_clock();

// Returns mpz_class random number of specified length
mpz_class mpz_rand(size_t length);

// Counts the number of lines in a given file
unsigned count_lines(std::string filename);

// Returns true if the specified file exists
bool file_exists(std::string filename);

// Estimate peak RAM usage from specified fraction size (given in MBs)
size_t peak_ram_usage(size_t fraction_size);

// Get recommended fraction split size based on available RAM
size_t recommended_split_size();


// Converts shorthand literals to a number (e.g. 25k -> 25 000)
uint64_t parse_shorthand_num(std::string str);

// Returns a string containing information regarding supported shorthand literals
std::string print_shorthand_info();

// Expands the given time (in ms) to days, hours, minutes, seconds, ms...
std::string format_time(double ms);

// Returns the total RAM in a system in MBs
size_t getTotalSystemMemoryMB();

#endif
