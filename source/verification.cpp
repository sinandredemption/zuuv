#include "verification.h"

#include "disk_mpz.h"
#include "continuant.h"
#include <string>
#include <sstream>
#include <filesystem>
void calc_continuant_from_file(std::string filename) {
	std::ifstream cf_terms_file(filename);
	CFTerms cf_terms;

	if (!cf_terms_file.good()) {
		throw "Can't open file '" + filename + "'";
		return;
	}

	std::stringstream buffer;
	buffer << cf_terms_file.rdbuf();

	uint64_t term;
	while (true)
	{
		if (!(buffer >> term)) break;
		else cf_terms.push_back(term);
	}

	mpz_class p_k = continuant(0, cf_terms.size(), cf_terms);
	// Offload p_k
//	disk_mpz disk_p_k(filename + "_p_k", disk_mpz::SplitSize);

}