#include "cruncher.h"
#include "benchmark.h"
#include <iostream>
#include <thread>
#include <filesystem>

int main() {
	std::cout << "Zuuv - Multi-precision Floating-point to Continued Fraction Cruncher" << std::endl;
	std::cout << "Version: v0.0 ALPHA (released 4/Mar/21) | ";

	std::cout << std::thread::hardware_concurrency() / 2
		<< " cores | hyperthreading disabled" << std::endl;

	std::filesystem::create_directory("iterations");
	std::filesystem::create_directory("disk_mpz");

	prompt:
	std::cout << std::endl;
	std::cout << "1 - Start a new computation" << std::endl;
	std::cout << "2 - Resume an existing computation" << std::endl;
	std::cout << "3 - Benchmark Continued Fraction Cruncher" << std::endl;
	std::cout << "4 - Benchmark Continuant Cruncher" << std::endl;
	std::cout << "\nEnter your choice: ";
	int choice;
	std::cin >> choice;

	try {
		switch (choice)
		{
		case 1: case 2:
		{
			std::string file = "";
			if (choice == 1) {
				std::cout << "\nEnter filename: ";
				std::cin >> file;
			}

			std::cout << "Enter # terms: ";
			uint64_t terms = 0;
			std::cin >> terms;

			std::cout << "Enter bytes per file (enter 0 for an all-in-RAM computation): ";
			size_t bytes_per_file = 0;
			std::cin >> bytes_per_file;

			if (bytes_per_file > 0) {
				//std::cout << "Enter number of threads to use: ";
				size_t nthreads = std::thread::hardware_concurrency() / 2;
				//std::cin >> nthreads;

				crunch_reg_cf_terms_on_disk(file, terms, bytes_per_file, nthreads);
			}
			else crunch_reg_cf_terms(file, terms);

			break;
		}

		case 3:
		{
			bool verify = false;
			do {
				std::cout << "Do you also want to verify the computations (yes/no)? ";
				std::string response = "";
				std::cin >> response;
				if (response == "yes") { verify = true; break; }
				if (response == "no") { verify = false; break; }
			} while (true);

			benchmark(verify);

			break;
		}
		case 4:
		{
			benchmark_continuant();
		}
		default:
			goto prompt;
			break;
		}
	}
	catch (const char* ex) {
		std::cerr << "\nERROR: " << *ex << "\nTerminating program..." << std::endl;
	}

	return 0;
}