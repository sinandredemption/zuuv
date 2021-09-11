#ifndef INC_MUL_H
#define INC_MUL_H
#include <ymp/ymp.h>
#include <mpir.h>
namespace multiplication {
  void init();

  size_t allocate_mem(size_t mulsize, size_t threads = 0);
  void release_mem();

  void mul(mpz_t res, const mpz_t s1, const mpz_t s2, size_t threads = 0);
  void addmul(mpz_t res, const mpz_t s1, const mpz_t s2, size_t threads = 0);
  void submul(mpz_t res, const mpz_t s1, const mpz_t s2, size_t threads = 0);

  extern size_t scratch_mem_size;
}

#endif  // INC_MUL_H