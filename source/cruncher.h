#ifndef INC_CONVERGENT_H
#define INC_CONVERGENT_H
#include <vector>
#include <mpir.h>
#include <mpirxx.h>
#include <string>

static_assert(sizeof(mp_limb_t) == sizeof(uint64_t), "Zuuv is currently compaitible only with 8-byte limbs");

typedef std::vector<mp_limb_t> CFTerms;

struct CFTermList
{
	bool updated = false;
	CFTerms list;

	size_t terms_on_disk = 0, idx = 0;

	// p_k  = continuant(list[0...k])
	// p_k1 = continuant(list[0...k-1])
	// q_k  = continuant(list[1...k])
	// q_k1 = continuant(list[1...k-1])
	struct Continuants
	{
		mpz_class p_k, p_k1, q_k, q_k1;
	} continuants;

	// updates peripheral constants for the list
	void update_continuants(bool single_term = false);

	void offload();	// Write to disk and free RAM
	inline std::string get_filename() const { return std::string("cf_terms") + std::to_string(idx); }
	void onload();
	void append(const CFTermList& cfterms);
};

CFTerms    basecase_calc_reg_cf_expansion(const mpz_class& num, const mpz_class& den, bool half);
CFTermList reg_cf_expansion(const mpz_class& num, const mpz_class& den, bool calc_convergents = false, bool half = false);

int perform_corrections(mpz_class& num, mpz_class& den, CFTermList& c);

void crunch_reg_cf_expansion(std::string file, size_t terms);
void crunch_reg_cf_expansion_on_disk(std::string file, size_t terms, size_t nthreads = 1);

#endif
