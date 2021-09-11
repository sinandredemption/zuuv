#include "mul.h"
#include <thread>
#include <cassert>
#include <mpirxx.h>

namespace multiplication
{
  void* scratch_mem_ptr;
  size_t scratch_mem_size;

  void init()
  {
    ymp::ensure_global_table<mp_limb_t>();
    scratch_mem_size = 0;
  }

  size_t allocate_mem(size_t mulsize, size_t threads)
  {
    if (threads == 0)
      threads = std::thread::hardware_concurrency();

    size_t memsize = ymp::LowLevel::mul_iPsize<mp_limb_t>(mulsize, threads);

    if (memsize > scratch_mem_size) {
      if (scratch_mem_size != 0)
        ymp::AlignedFree(scratch_mem_ptr);
      scratch_mem_ptr = ymp::AlignedMalloc(memsize);

      if (!scratch_mem_ptr)
        throw "Can't allocate enough memory";

      scratch_mem_size = memsize;
    }

    return memsize;
  }

  void release_mem()
  {
    ymp::AlignedFree(scratch_mem_ptr);
    scratch_mem_size = 0;
  }

  void mul(mpz_t res, const mpz_t s1, const mpz_t s2, size_t threads)
  {
#ifndef FORCE_MPIR_MUL
    if (threads == 0)
      threads = std::thread::hardware_concurrency();

    if (s1->_mp_size == 0 || s2->_mp_size == 0)
    {
      mpz_set_ui(res, 0);
      return;
    }

    assert((s1->_mp_size) > 0 && (s2->_mp_size > 0));

    size_t mulsize = mpz_size(s1) + mpz_size(s2);

    allocate_mem(mulsize, threads);

    if (res->_mp_alloc < mulsize) {
      mpz_realloc(res, mulsize);
    }

    std::memset(res->_mp_d, 0, res->_mp_alloc * sizeof(mp_limb_t));

    ymp::BasicParameters params(ymp::get_global_table(), scratch_mem_ptr, scratch_mem_size, threads);

    ymp::LowLevel::mul(params, res->_mp_d, s1->_mp_d, s1->_mp_size, s2->_mp_d, s2->_mp_size);

    res->_mp_size = mulsize;
    if (res->_mp_d[mulsize - 1] == 0)
      res->_mp_size--;
#else
    mpz_mul(res, s1, s2);
#endif
  }

  void addmul(mpz_t res, const mpz_t s1, const mpz_t s2, size_t threads)
  {
    mpz_class tmp;
    mul(tmp.get_mpz_t(), s1, s2, threads);
    mpz_add(res, res, tmp.get_mpz_t());
  }

  void submul(mpz_t res, const mpz_t s1, const mpz_t s2, size_t threads)
  {
    mpz_class tmp;
    mul(tmp.get_mpz_t(), s1, s2, threads);
    mpz_sub(res, res, tmp.get_mpz_t());
  }
}