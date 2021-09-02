#include "cruncher.h"
#include "benchmark.h"
#include <iostream>
#include <thread>
#include <fstream>
#include <filesystem>
#include <sstream>

int main() {
  std::cout << "Zuuv [EXPERIMENTAL] - Multi-precision Floating-point to Continued Fraction Cruncher" << std::endl;
  std::cout << "Version: v1.0 BETA (released 11/Apr/21) | Syed Fahad ( sydfhd AT gmail.com )\n";

  std::cout << "Using " << std::thread::hardware_concurrency() << " threads" << std::endl;

  std::filesystem::create_directory("iterations");
  std::filesystem::create_directory("disk_mpz");

prompt:
  std::cout << std::endl;

  std::cout << "1 - Start a new computation"                                 << std::endl;
  std::cout << "2 - Resume an existing computation"             << std::endl << std::endl;

  std::cout << "3 - Benchmark continued fraction cruncher"                   << std::endl;
  std::cout << "4 - Benchmark continuant cruncher"                           << std::endl;
  std::cout << "5 - Benchmark disk-based multiplication"                     << std::endl;
  std::cout << "6 - Benchmark RAM-based multiplication"         << std::endl << std::endl;

  std::cout << "7 - View incrementally largest terms in computed iterations" << std::endl;

  std::cout << "\nEnter your choice: ";
  int choice;
  std::cin >> choice;

  try {
    switch (choice)
    {
    case 1:
    {
      size_t terms = 0;
      while (terms == 0) {
        std::cout << "Enter # terms: ";
        std::cin >> terms;
      }

      int ram_based = 0;

      std::cout << "1 - (Default) RAM based computation (~" << (terms * 1.668041967e-5) << "MB required)\n";
      std::cout << "2 - Disk based computation" << std::endl;
      std::cout << "Choose: ";
      std::cin >> ram_based;

      size_t bytes_per_file = 0;
      if (ram_based == 2) {
        double file_size = 0;

        std::cout << "\nPeak RAM usage ~= 45x file size" << std::endl;
        std::cout << "Enter file size (MBs, default = 10): ";
        std::cin >> file_size;
        bytes_per_file = file_size * 1024. * 1024. + 1;
      }

      std::string file = "";
      if (choice == 1) {
        std::cout << "\nEnter filename: ";
        std::cin >> file;
      }

      if (bytes_per_file > 0) {
        std::ofstream options("zuuv.txt");
        options << bytes_per_file << std::endl;
        options << terms << std::endl;
        options.close();

        std::cout << "WARNING: The file must contain hexadecimal digits.\n"
          << "Decimal digits for disk-based computations are not yet supported.\n" << std::endl;
        size_t nthreads = std::thread::hardware_concurrency();

        crunch_reg_cf_terms_on_disk(file, terms, bytes_per_file, nthreads);
      }
      else crunch_reg_cf_terms(file, terms);

      break;
    }

    case 2:
    {
      size_t bytes_per_file, terms, nthreads = std::thread::hardware_concurrency() / 2;

      std::ifstream options("zuuv.txt");
      if (!options.good()) throw "Can't load options from file 'zuuv.txt'";
      options >> bytes_per_file >> terms;
      options.close();

      crunch_reg_cf_terms_on_disk("", terms, bytes_per_file, nthreads);

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

      benchmark_cf_cruncher(verify);

      break;
    }
    case 4:
    {
      benchmark_continuant();
      break;
    }
    case 5:
    {
      benchmark_disk_mul();
      break;
    }
    case 6:
    {
      benchmark_ram_mul();
      break;
    }
    case 7:
    {
      std::cout << std::endl;
      std::cout << "\tTerm | Position" << std::endl;
      unsigned long long int max = 0, term_n = 1;
      
      for (int i = 1;; ++i) {
        if (!std::filesystem::exists("iterations/iteration" + std::to_string(i) + ".txt"))
          break;

        std::fstream file("iterations/iteration" + std::to_string(i) + ".txt");
        std::stringstream buffer;
        buffer << file.rdbuf();

        for (unsigned long long int term; buffer >> term; term_n++)
          if (term > max) {
            max = term;
            std::cout << "\t" << max << " | " << term_n << std::endl;
          }
      }

      std::cout << std::endl;
      break;
    }
    default:
      goto prompt;
      break;
    }
  }
  catch (const char* ex) {
    std::cerr << "\nERROR: " << ex << "\nTerminating program..." << std::endl;
  }

  std::cout << "Press enter to continue...";
  std::cin.get(); std::cin.get();
  return 0;
}