#include "utils.h"
#include "params.h"
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
	for (int i = 0; i < length; ++i)
	{
		static uint64_t x = 0x314159265358979ULL;

		x ^= x << 13;
		x ^= x >> 7;
		x ^= x << 17;

		ptr[i] = x;
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

size_t peak_ram_usage(size_t split_size)
{
	return split_size * Params::RAMUsagePerMBofFraction + Params::BaselineRAMUsage;
}

size_t recommended_split_size()
{
	return size_t(0.9 * double(getTotalSystemMemoryMB() - Params::BaselineRAMUsage)
		                            / Params::RAMUsagePerMBofFraction);
}

uint64_t parse_shorthand_num(std::string str)
{
	std::stringstream ss;
	ss << str.substr(0, str.size() - 1);

	if (std::isalpha(str[str.size() - 1]))
	{
		double n;
		ss >> n;

		switch (str[str.size() - 1])
		{
		case 'K':
		case 'k':
			return static_cast<uint64_t>(n * 1e3);
			break;
		case 'm':
		case 'M':
			return static_cast<uint64_t>(n * 1e6);
			break;
		case 'b':
		case 'B':
			return static_cast<uint64_t>(n * 1e9);
			break;
		case 't':
		case 'T':
				return static_cast<uint64_t>(n * 1e12);
				break;
		default:
			return 0;
			break;
		}
	}
	else
	{
		uint64_t n;
		ss >> n;

		return n;
	}
}

std::string print_shorthand_info()
{
	std::stringstream ss;

	ss << "k = K = 1 000             (10^3)" << std::endl;
	ss << "m = M = 1 000 000         (10^6)" << std::endl;
	ss << "b = B = 1 000 000 000     (10^9)" << std::endl;
	ss << "t = T = 1 000 000 000 000 (10^12)" << std::endl;

	return ss.str();
}

std::string format_time(double time_ms)
{
	std::stringstream ss;
	int d, h, m, s, ms;
	d = h = m = s = ms = 0;

	if (time_ms < 1000)
	{
		ms = (int)time_ms;
	}
	else if (time_ms < 60 * 1000)	// One minute
	{
		s = time_ms / 1000;
		ms = time_ms - 1000 * s;
	}
	else if (time_ms < 60 * 60 * 1000)	// One hour
	{
		s = time_ms / 1000;
		ms = time_ms - 1000 * s;

		m = s / 60;
		s -= m * 60;
	}
	else if (time_ms < 24 * 60 * 60 * 1000) // One day
	{
		s = time_ms / 1000;
		ms = time_ms - 1000 * s;

		int m = s / 60;
		s -= m * 60;

		h = m / 60;
		m -= h * 60;
	}
	else
	{
		s = time_ms / 1000;
		ms = time_ms - 1000 * s;

		int m = s / 60;
		s -= m * 60;

		h = m / 60;
		m -= h * 60;

		d = h / 24;
		h -= 24 * d;
	}

	if (d)  ss << d  << "d ";
	if (h)  ss << h  << "h ";
	if (m)  ss << m  << "m ";
	if (s)  ss << s  << "s ";
	if (ms) ss << ms << "ms";

	return ss.str();
}

// Based on shameless copy paste from stackoverflow
#ifdef _WIN32
#include <windows.h>

unsigned long long getTotalSystemMemoryMB()
{
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return status.ullTotalPhys / (1024*1024);
}
#elif defined(linux)
#include <unistd.h>

unsigned long long getTotalSystemMemoryMB()
{
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return (pages * page_size) / (1024*1024);
}
#endif
