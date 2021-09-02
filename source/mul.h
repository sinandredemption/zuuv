#ifndef INC_MUL_H
#define INC_MUL_H
#include <ymp/ymp.h>
#include <mpir.h>
namespace multiplication {
  void init();
  void allocate(size_t bytes);
  void mul(mpz_t res, const mpz_t s1, const mpz_t s2, size_t threads);
}

#endif  // INC_MUL_H