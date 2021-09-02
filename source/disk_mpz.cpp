#include "disk_mpz.h"
#include "utils.h"
#include "mul.h"
#include <fstream>
#include <iostream>
#include <exception>
#include <thread>
#include <filesystem>

/*
  The default constructor
  Params: "name" = name of the fraction (used in disk i/o), "filename" = path to file containing
  hex digits, "bytes_per_file" = terms_on_disk of each file
  Note that the constructor ignores the integer part of fraction, i.e. 3.14159... will be
  internally stored as 0.14159265...
*/
disk_mpq::disk_mpq(std::string name, std::string filename, size_t bytes_per_file)
  : num(name + "_num", bytes_per_file), den(name + "_den", bytes_per_file)
{
  if (filename == "") { // Find files on disk
    size_t n = 0;
    do {
      std::string f_num("disk_mpz\\" + name + "_num_" + std::to_string(n)),
        f_den("disk_mpz\\" + name + "_den_" + std::to_string(n));
      if (file_exists(f_num))
        num.pushback_file(name + "_num_" + std::to_string(n));
      if (file_exists(f_den))
        den.pushback_file(name + "_den_" + std::to_string(n));
      if (!file_exists(f_num) && !file_exists(f_den)) break;
      ++n;
    } while (true);
    return;
  }

  std::cerr << "Reading from file: " << filename << "...";
  try {
    std::ifstream inp(filename);
    if (!inp) throw "Couldn't open file";

    // Jump to end of file... Read in reverse while < bytes_per_file digits left
    inp.seekg(0, std::ios::end);

    char* buffer = new char[bytes_per_file * 2 + 1]; // Buffer to store the digits read

    for (size_t n = 0; (size_t) inp.tellg() > bytes_per_file * 2 + 1; ++n) { // "bytes_per_file + 1" to ignore the dot
      inp.seekg(-(int64_t)(bytes_per_file*2), std::ios::cur);

      // Read digits into buffer
      inp.read(buffer, (bytes_per_file * 2));
      buffer[bytes_per_file * 2] = '\0';

      // Pushback the read digits
      num.pushback_digits(buffer);
      den.pushback_zero();

      inp.seekg(-(int64_t)(bytes_per_file*2), std::ios::cur);
    }

    // Read the remaining digits (the first few digits)
    if (inp.tellg() > 2) {  // Ignore the decimal point and anything before it
      std::memset(buffer, 0, sizeof(char) * ((2*bytes_per_file) + 1));

      size_t offset = inp.tellg();
      inp.seekg(2, std::ios::beg);  // Ignore the decimal point...
      inp.read(buffer, offset - 2);

      // Pushback the read digits
      num.pushback_digits(buffer);

      // Calculate denominator
      mpz_class denom(1);
      denom <<= 4 * (offset - 2); // 4 bits per character
      den.pushback_mpz(denom);
    }

    delete[] buffer;
  }
  catch (const char* msg) {
    std::cerr << "FAILED: " << msg << std::endl;
  }
  catch (std::string msg) {
    std::cerr << "FAILED: " << msg << std::endl;
  }

  std::cout << "LOADED" << std::endl;
}

void disk_mpz::canonicalize()
{
  // Truncate leading zereos
  while (get_top_mpz() == mpz_class(0))
    pop_top();

  // Process carries
  const mpz_class Base(mpz_class(1) << (8 * bytes_per_file));
  if (get_top_mpz() < 0) {
    set_neg();
    bool carry = false;
    for (int j = 0; j < files(); ++j) {
      mpz_class m_j (-mpz_at(j));
      if (carry) m_j -= mpz_class(1);

      carry = (m_j < 0);
      if (carry) m_j += Base;

      assert(m_j >= 0);
      set_mpz(j, m_j);
    }
    assert(carry == 0);
  }
  else {
    bool carry = false;
    for (int j = 0; j < files(); ++j) {
      mpz_class m_j = mpz_at(j);
      if (carry) m_j -= 1;
      carry = (m_j < 0);
      if (carry) m_j += Base;

      set_mpz(j, m_j);
    }

    assert(carry == 0);
  }

  // Truncate leading zereos
  while (get_top_mpz() == mpz_class(0))
    pop_top();
}

void disk_mpz::mt_canonicalize()
{
  // Update file list
  digit_files.clear();
  for (size_t j = 0; file_exists(get_filename(j)); ++j)
    digit_files.push_back(get_filename(j));

  // Truncate leading zereos
  while (get_top_mpz() == mpz_class(0))
    pop_top();

  // Process carries
  mpz_class carry = 0;

  for (int j = 0; j < files(); ++j) {
    mpz_class m_j(mpz_at(j));
    m_j += carry;

    mpz_tdiv_q_2exp(carry.get_mpz_t(), m_j.get_mpz_t(), 8 * bytes_per_file);
    mpz_tdiv_r_2exp(m_j.get_mpz_t(), m_j.get_mpz_t(), 8 * bytes_per_file);

    set_mpz(j, m_j);
  }

  if (carry != 0)
    pushback_mpz(carry);

  canonicalize();
}

// Multiply by m, and return result
disk_mpz disk_mpz::mul(mpz_class m, std::string mul_name) const
{
  disk_mpz result(mul_name, bytes_per_file);

  mpz_class carry(0);
  for (auto digit_file : digit_files) {
    assert(carry < m);

    mpz_class m2 = read_mpz(digit_file) * m + carry;
    // carry = m2 / 2 ^ (8 * bytes_per_file)
    mpz_tdiv_q_2exp(carry.get_mpz_t(), m2.get_mpz_t(), 8 * bytes_per_file);
    // m2 % 2 ^ (8 * bytes_per_file)
    mpz_tdiv_r_2exp(m2.get_mpz_t(), m2.get_mpz_t(), 8 * bytes_per_file);

    result.pushback_mpz(m2);
  }

  if (carry > 0)
    result.pushback_mpz(carry);

  return result;
}


disk_mpz disk_mpz::cross_mult_sub(std::string name,
  const disk_mpz& a, const mpz_class& a1, const disk_mpz& b, const mpz_class& b1,
  size_t bytes_per_file, size_t threads)
{
  disk_mpz res(name, bytes_per_file);

  mpz_class carry(0), m, tmp;
  for (size_t i = 0; i < std::max(a.files(), b.files()); ++i) {
    // m = (a.mpz_at(i) * a1 - b.mpz_at(i) * b1) + carry;
    m = carry;
    //mpz_addmul(m.get_mpz_t(), a.mpz_at(i).get_mpz_t(), a1.get_mpz_t());
    multiplication::mul(tmp.get_mpz_t(), a.mpz_at(i).get_mpz_t(), a1.get_mpz_t(), threads);
    m += tmp;
    //mpz_submul(m.get_mpz_t(), b.mpz_at(i).get_mpz_t(), b1.get_mpz_t());
    multiplication::mul(tmp.get_mpz_t(), b.mpz_at(i).get_mpz_t(), b1.get_mpz_t(), threads);
    m -= tmp;

    // carry = m / 2 ^ (8 * bytes_per_file)
    //     m = m % 2 ^ (8 * bytes_per_file)
    mpz_tdiv_q_2exp(carry.get_mpz_t(), m.get_mpz_t(), 8 * bytes_per_file);
    mpz_tdiv_r_2exp(    m.get_mpz_t(), m.get_mpz_t(), 8 * bytes_per_file);

    res.pushback_mpz(m);
  }

  if (carry != 0)
    res.pushback_mpz(carry);

  res.canonicalize();

  return res;
}

disk_mpz disk_mpz::sub(const disk_mpz& s, std::string sub_name) const
{
  disk_mpz result(sub_name, bytes_per_file);


  // Find the greater among two
  const disk_mpz *greater = NULL, *lesser = NULL;
  for (int i = (int)std::max(files(), s.files()) - 1; i >= 0; --i) {
    int cmp = mpz_cmp(mpz_at(i).get_mpz_t(), s.mpz_at(i).get_mpz_t());
    if (cmp == 0) continue;
    if (cmp < 0) {
      greater = &s, lesser = this;
      break;
    }
    if (cmp > 0) {
      greater = this, lesser = &s;
      break;
    }
  }

  // Perform subtraction
  const mpz_class Base(mpz_class(1) << (8 * bytes_per_file));
  bool carry = false;
  for (size_t j = 0; j < std::max(files(), s.files()); ++j) {
    mpz_class a = greater->mpz_at(j) - lesser->mpz_at(j);
    if (carry) a -= 1;
    carry = a < 0;
    if (carry) a += Base;

    result.pushback_mpz(a);
  }

  assert(carry == 0);
  if (greater != this)
    result.set_neg();

  // Truncate leading zeroes
  while (result.get_top_mpz() == mpz_class(0))
    result.pop_top();

  return result;
}

void disk_mpz::pushback_mpz(const mpz_class& mpz)
{
  set_mpz(digit_files.size(), mpz);
}

void disk_mpz::set_mpz(size_t j, const mpz_class& mpz, bool update_file_list)
{
  assert(mpz_size(mpz.get_mpz_t()) < bytes_per_file * sizeof(mp_limb_t));

  std::string filename;
  if (j < digit_files.size()) {
    filename = digit_files[j];
    remove(filename.c_str());
  }
  else {
    assert(j == digit_files.size() || !update_file_list);
    filename = get_filename(j);
    if (update_file_list) // In case of multi-threaded access, don't update file list
      digit_files.push_back(filename);
  }

  std::ofstream f(filename, std::ios::out | std::ios::binary);
  if (!f.good())
    throw name + ".set_mpz(): Can't open file '" + filename + "'";

  f.write((const char*)&(mpz.get_mpz_t()->_mp_size), sizeof(mpz.get_mpz_t()->_mp_size));
  f.write((const char*)(mpz.get_mpz_t()->_mp_d),
    sizeof(*mpz.get_mpz_t()->_mp_d) * std::abs(mpz.get_mpz_t()->_mp_size));
  f.close();
}

void disk_mpz::pushback_digits(char* buffer)
{
  mpz_t temp_mpz;
  mpz_init(temp_mpz);

  if (mpz_set_str(temp_mpz, buffer, 16) != 0) {
    throw name + ".pushback_digits(): set_str() failed";
  }

  std::string filename = get_filename(digit_files.size());
  digit_files.push_back(filename);

  std::ofstream f(filename, std::ios::out | std::ios::binary);
  if (!f.good())
    throw name + ".pushback_digits(): Can't open file '" + filename + "'";

  f.write((char*)&(temp_mpz->_mp_size), sizeof(temp_mpz->_mp_size));
  f.write((char*)(temp_mpz->_mp_d), sizeof(*temp_mpz->_mp_d) * std::abs(temp_mpz->_mp_size));
  f.close();

  mpz_clear(temp_mpz);
}

void disk_mpz::pushback_zero()
{
  char str[] = "0";
  pushback_digits(str);
}

void disk_mpz::pushback_file(std::string filename)
{
  digit_files.push_back("disk_mpz\\" + filename);
}

mpz_class disk_mpz::get_top2_mpz() const
{
  if (files() == 1) return get_top_mpz();
  return (get_top_mpz() << (8 * bytes_per_file)) | mpz_at(digit_files.size() - 2);
}
mpz_class disk_mpz::to_mpz() const
{
  mpz_class out(0);
  for (size_t i = 0; i < files(); ++i)
    out += mpz_at(i) << (8 * bytes_per_file * i);
  if (sgn < 0) out = -out;
  return out;
}

mpz_class disk_mpz::value() const
{
  mpz_class out(0);
  for (size_t i = 0; i < files(); ++i)
    out += mpz_at(i) << (8 * bytes_per_file * i);
  return out;
}

disk_mpz::~disk_mpz()
{
}

void disk_mpz_move(disk_mpz& to, disk_mpz& from)
{
  while (to.files() > 0)
    to.pop_top();

  for (auto f : from.digit_files) {
    std::filesystem::rename(f, to.get_filename(to.files()));
    to.digit_files.push_back(to.get_filename(to.files()));
  }

  to.sgn = from.sgn;
}

mpz_class disk_mpz::read_mpz(std::string filename)
{
  std::ifstream bin_in(filename, std::ios::in | std::ios::binary);

  mpz_class mpz;

  int newsize;
  assert(sizeof(mpz.get_mpz_t()->_mp_size) == sizeof(newsize));
  bin_in.read((char*)&newsize, sizeof(newsize));
  mpz_realloc(mpz.get_mpz_t(), (mp_size_t)std::abs(newsize));

  bin_in.read((char*)(mpz.get_mpz_t()->_mp_d), sizeof(*mpz.get_mpz_t()->_mp_d) * std::abs(newsize));
  mpz.get_mpz_t()->_mp_size = newsize;

  bin_in.close();

  return mpz;
}
