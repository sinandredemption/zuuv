#ifndef INC_DISK_MPZ_H
#define INC_DISK_MPZ_H
#include <mpir.h>
#include <mpirxx.h>
#include <string>
#include <vector>
#include <cassert>
#include <fstream>
#include <mutex>

/*
  An disk_mpz is represented by a vector of files (vector<string> digit_files, where the string
  holds the path to file containing an mpz dump of the digits), with "bytes_per_file" number of
  digits per file. It evaluates to:
  digit_file[0] * 2^(bytes_per_file*0) + digit_file[1] * 2^(bytes_per_file*1) + ...
*/
class disk_mpz
{
  std::vector<std::string> digit_files;
  std::string name;
  size_t bytes_per_file;
  int sgn;

public:
  const char* SubDir = "disk_mpz/";

  void pushback_mpz(const mpz_class& mpz);
  void set_mpz(size_t j, const mpz_class& mpz, bool update_file_list = true);

  // Pushback digits from buffer continaing hex string of digits into a NEW file
  void pushback_digits(char *buffer);
  void pushback_zero(); // Pushback new file containing zero x bytes_per_file
  void pushback_file(std::string filename);

  void canonicalize();  // Truncate leading zereos and process boolean carries
  void mt_canonicalize(); // Canonicalize for Multi-threaded procedures with carries often > 1
  disk_mpz mul(mpz_class, std::string) const;
  disk_mpz sub(const disk_mpz&, std::string) const;

  // Compute a * a1 - b * b1
  static disk_mpz cross_mul_sub(std::string name,
    const disk_mpz& a, const mpz_class& a1, const disk_mpz& b, const mpz_class& b1,
    size_t bytes_per_file, size_t threads = 0);

  std::string getname() const { return name; }
  std::string get_filename(size_t n) const { return std::string(SubDir) + name + "_" + std::to_string(n); }
  size_t files() const { return digit_files.size(); }
  mpz_class mpz_at(size_t n) const { return n >= digit_files.size() ? mpz_class(0) : read_mpz(digit_files[n]); }
  void set_neg() { sgn = -1; }
  void set_pos() { sgn = 1; }
  int sign() const { return sgn; }
  std::string get_top() const { return digit_files[digit_files.size() - 1]; }
  mpz_class get_top_mpz() const { return read_mpz(digit_files.back()); }
  mpz_class get_top2_mpz() const;
  void pop_top() { remove(get_top().c_str()); digit_files.pop_back(); }
  mpz_class to_mpz() const; // Useful for debugging
  mpz_class value() const;  // Useful for debugging

  void destroy() { while (digit_files.size()) pop_top(); }

  disk_mpz(std::string name_, size_t bytes_per_file_)
    : name(name_), bytes_per_file(bytes_per_file_), sgn(0) {}
  ~disk_mpz();

  // --- Helper functions ---
  friend void disk_mpz_move(disk_mpz&, disk_mpz&);
  static mpz_class read_mpz(std::string filename);
};

/*
  An disk_mpq class is represented by a disk_mpz numerator and denomenator which share a common prefix
*/
class disk_mpq
{
  disk_mpz num, den;
public:
  disk_mpq(std::string, std::string, size_t); // Name, file, digits per file

  disk_mpz& get_num() { return num; }
  const disk_mpz& get_num() const { return num; }
  disk_mpz& get_den() { return den; }
  const disk_mpz& get_den() const { return den; }

  void destroy() { num.destroy(); den.destroy(); }
  // -- DEBUG
  double get_d() const {
    return mpq_class(disk_mpz::read_mpz(num.get_top()), disk_mpz::read_mpz(den.get_top())).get_d();
  }
  // --
};

#endif
