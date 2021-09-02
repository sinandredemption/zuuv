#include "mul.h"

namespace multiplication {
  void* scratch_mem_ptr;
  size_t scratch_mem_size;

  void init() {
    ymp::ensure_global_table<mp_limb_t>();
    scratch_mem_size = 0;
  }

  void allocate_mem(size_t bytes) {
    if (bytes > scratch_mem_size) {
      if (scratch_mem_size != 0)
        ymp::AlignedFree(scratch_mem_ptr);

      scratch_mem_ptr = ymp::AlignedMalloc(bytes);
      scratch_mem_size = bytes;
    }
  }

  void release_mem() {
    ymp::AlignedFree(scratch_mem_ptr);
    scratch_mem_size = 0;
  }

  void mul(mpz_t res, const mpz_t s1, const mpz_t s2, size_t threads) {
    size_t mulsize = mpz_size(s1) + mpz_size(s2);
    size_t memsize = ymp::LowLevel::mul_iPsize<mp_limb_t>(mulsize, threads);

    allocate_mem(memsize);

    if (res->_mp_alloc < mulsize) {
      mpz_realloc(res, mulsize);
      std::memset(res->_mp_d, 0, mulsize * sizeof(mp_limb_t));
    }

    ymp::BasicParameters params(ymp::get_global_table(), scratch_mem_ptr, scratch_mem_size, threads);

    ymp::LowLevel::mul(params, res->_mp_d, s1->_mp_d, s1->_mp_size, s2->_mp_d, s2->_mp_size);

    res->_mp_size = mulsize;
  }
}