#include <ymp/ymp.h>
#include <mpir.h>
#include <iostream>
#include <thread>
#include <ymp/LowLevel.h>
void mpir_mul(uint64_t *a, uint64_t *b, uint64_t *c, size_t N) {
  mpn_mul_n(c, a, b, N);
}

void ymp_mul(uint64_t* a, uint64_t* b, uint64_t* c, size_t N) {
  static auto lookup = ymp::MakeLookupTable().get();

  size_t scratch_mem_size = ymp::LowLevel::mul_iPsize<uint64_t>(N, std::thread::hardware_concurrency());

  auto uptr = ymp::SmartPointer<>::malloc_uptr(scratch_mem_size, ymp::DEFAULT_ALIGNMENT);

  void* scratch_mem = uptr.get();

  ymp::BasicParameters params(lookup, scratch_mem, scratch_mem_size, std::thread::hardware_concurrency());

  ymp::LowLevel::mul(params, c, a, N, b, N);
}

int main() {
  size_t mulsize = 0;

  std::cout << "mulsize: ";
  std::cin >> mulsize;

  uint64_t* a, * b, * c1, * c2;
  a = new uint64_t[mulsize];
  b = new uint64_t[mulsize];
  c1 = new uint64_t[mulsize * 2];
  c2 = new uint64_t[mulsize * 2];

  for (uint64_t* a_ = a; a_ != a + mulsize; ++a_) {
    *a_ = uint64_t(rand()) * uint64_t(rand())* uint64_t(rand())* uint64_t(rand())* uint64_t(rand());
  }


  for (uint64_t* b_ = b; b_ != b + mulsize; ++b_) {
    *b_ = uint64_t(rand()) * uint64_t(rand()) * uint64_t(rand()) * uint64_t(rand()) * uint64_t(rand());
  }

  mpir_mul(a, b, c1, mulsize);
  ymp_mul(a, b, c2, mulsize);

  std::cin.get();
  return 0;
}