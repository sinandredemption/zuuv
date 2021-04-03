#ifndef INC_CONTINUANT_H
#define INC_CONTINUANT_H
#include "cruncher.h"
#include "params.h"
#include <mpir.h>
#include <mpirxx.h>
#include <map>
#include <cassert>

typedef struct { mpz_class first, second; } mpz_pair;

namespace ContinuantCache {
	mpz_class cached_continuant(size_t s, size_t t, size_t mid, const CFTerms&);
	mpz_class parallel_cached_continuant(size_t s, size_t t, size_t mid, const CFTerms&);
	void clear();
}

mpz_class continuant(size_t s, size_t t, const CFTerms&, size_t split_point = 0);
mpz_class parallel_continuant(size_t s, size_t t, const CFTerms&, size_t split_point = 0);

// TODO Replace mpz_addmul with mpn_addmul
// Returns { continuant(s, t), continuant(s, t - 1) }
inline mpz_pair continuant_pair_right(size_t s, size_t t, const CFTerms& terms) {
	mpz_class K1, K2;
	auto k_n = &K1, k_n1 = &K2;

	mpz_realloc(k_n->get_mpz_t(), Params::DefaultContAllocSize);
	mpz_realloc(k_n1->get_mpz_t(), Params::DefaultContAllocSize);

	K1 = terms[s++], K2 = 1;

	while (s != t) {
		std::swap(k_n, k_n1);
		mpz_addmul_ui(k_n->get_mpz_t(), k_n1->get_mpz_t(), terms[s++]);
	}

	mpz_addmul_ui(k_n1->get_mpz_t(), k_n->get_mpz_t(), terms[s]);
	return { *k_n1, *k_n };
}

// Returns { continuant(s, t), continuant(s + 1, t) }
inline mpz_pair continuant_pair_left(size_t s, size_t t, const CFTerms& terms) {
	mpz_class K1, K2;
	auto k_n = &K1, k_n1 = &K2;

	mpz_realloc(k_n->get_mpz_t(), Params::DefaultContAllocSize);
	mpz_realloc(k_n1->get_mpz_t(), Params::DefaultContAllocSize);

	K1 = terms[t--], K2 = 1;

	while (t != s) {
		std::swap(k_n, k_n1);
		mpz_addmul_ui(k_n->get_mpz_t(), k_n1->get_mpz_t(), terms[t--]);
	}

	mpz_addmul_ui(k_n1->get_mpz_t(), k_n->get_mpz_t(), terms[t]);
	return { *k_n1, *k_n };
}

inline mpz_class basic_continuant(size_t s, size_t t, const CFTerms& terms) {
	auto mid = (s + t) / 2;
	
	mpz_pair lo(continuant_pair_right(s, mid, terms));
	mpz_pair hi(continuant_pair_left(mid + 1, t, terms));

	return lo.first * hi.first + lo.second * hi.second;
}

#endif
