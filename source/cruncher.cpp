#include "cruncher.h"
#include "continuant.h"
#include "io.h"
#include "utils.h"
#include "disk_mpz.h"
#include "mul.h"
#include <cassert>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <thread>
#include <filesystem>

#define SLOPPY_BASECASE_DIVISION

CFTerms basecase_calc_reg_cf_expansion(const mpz_class& num, const mpz_class& den, bool half)
{
	CFTerms cfterms;
	cfterms.reserve(Params::DefaultCFTermsAlloc);

	mpz_class r, N(num), D(den), s;
	mpz_class* n = &N, * d = &D;

	if (half) s = sqrt(num);

	while (*d != 0)
	{
		auto longdivide_and_pushback = [&]() {
			mpz_tdiv_q(r.get_mpz_t(), n->get_mpz_t(), d->get_mpz_t());

			if (!r.fits_ui_p()) throw "Unusually large convergent";

			cfterms.push_back(r.get_ui());
			mpz_submul_ui(n->get_mpz_t(), d->get_mpz_t(), r.get_ui());
		};

#ifdef SLOPPY_BASECASE_DIVISION
		auto msl = [](const mpz_class& n)	// Most significant limb
		{
			return n.get_mpz_t()->_mp_d[mpz_size(n.get_mpz_t()) - 1];
		};

		if (mpz_size(d->get_mpz_t()) > 1 && mpz_size(n->get_mpz_t()) == mpz_size(d->get_mpz_t()))
		{
			auto n_top = msl(*n), d_top = msl(*d);

			auto r1 =      n_top  /      d_top;
			auto r2 = (1 + n_top) /      d_top;
			auto r3 =      n_top  / (1 + d_top);
			auto r4 = (1 + n_top) / (1 + d_top);

			if (r1 == r2 && r2 == r3 && r3 == r4)
			{
				cfterms.push_back(r1);
				mpz_submul_ui(n->get_mpz_t(), d->get_mpz_t(), r1);
			}
			else longdivide_and_pushback();
		}
		else longdivide_and_pushback();
#else
		longdivide_and_pushback();
#endif
		if (half && *n < s) break;

		std::swap(n, d);
	}

	if (cfterms.size() <= 1)
	{
		throw "Unusual behaviour: single-term continued fraction expansion";
	}

	return cfterms;
}

int perform_corrections(mpz_class& num, mpz_class& den, CFTermList& c)
{
	int corrections = 0;

	while (mpz_sgn(num.get_mpz_t()) < 0 || mpz_sgn(den.get_mpz_t()) < 0 || (num < den))
	{
		corrections++;
		mpz_swap(num.get_mpz_t(), den.get_mpz_t());

		if (mpz_sgn(den.get_mpz_t()) < 0)
		{
			assert(num > 0);

			mpz_neg(num.get_mpz_t(), num.get_mpz_t());
			mpz_neg(den.get_mpz_t(), den.get_mpz_t());
		}

		if (c.list.empty())
		{
			mp_limb_t lastterm;
			std::ifstream cfterms_file(c.get_filename(), std::ios::in | std::ios::binary);
			cfterms_file.seekg((c.terms_on_disk - corrections) * sizeof(mp_limb_t));
			cfterms_file.read((char*)&lastterm, sizeof(mp_limb_t));
			mpz_addmul_ui(num.get_mpz_t(), den.get_mpz_t(), lastterm);
		}
		else
		{
			mpz_addmul_ui(num.get_mpz_t(), den.get_mpz_t(), c.list.back());
			c.list.pop_back();

			if (c.list.empty()) return corrections;
		}
	}

	return corrections;
}

CFTermList reg_cf_expansion(const mpz_class& num, const mpz_class& den, bool calc_convergents, bool half)
{
	assert(num >= 0);
	assert(den > 0);

	size_t upper_half_size = std::min(mpz_size(num.get_mpz_t()), mpz_size(den.get_mpz_t())) / 2;

	if (upper_half_size < Params::ThresholdCFUseBasicProcedure)
	{
		CFTermList out;
		out.list = basecase_calc_reg_cf_expansion(num, den, half);

		if (calc_convergents) out.update_continuants();

		return out;
	}
	
	size_t shift = upper_half_size * sizeof(mp_limb_t) * 8;
	// Take CF of upper half
	CFTermList upperhalf_cf_terms(reg_cf_expansion(num >> shift, den >> shift, true, true));

	assert(upperhalf_cf_terms.terms_on_disk > 0 || !upperhalf_cf_terms.list.empty());
	
	if (upperhalf_cf_terms.list.size() > Params::CFTermsUseDisk)
		upperhalf_cf_terms.offload();
	
	// Determine correction fraction
	// Numerator = b * p_{k-1} - a * q_{k-1}
	// Denomenator = a * q_k - b * p_k
	assert(upperhalf_cf_terms.updated);

	// Numerator and denomenator of the correction factor
	mpz_class correction_num(0), correction_den(0);

  // correction_den = c_den1 - c_den2 (same for correction_num)
  mpz_class c_den1, c_den2, c_num1, c_num2;

  multiplication::mul(c_den1.get_mpz_t(), num.get_mpz_t(), upperhalf_cf_terms.continuants. q_k.get_mpz_t());
  multiplication::mul(c_den2.get_mpz_t(), den.get_mpz_t(), upperhalf_cf_terms.continuants. p_k.get_mpz_t());
  multiplication::mul(c_num1.get_mpz_t(), den.get_mpz_t(), upperhalf_cf_terms.continuants.p_k1.get_mpz_t());
  multiplication::mul(c_num2.get_mpz_t(), num.get_mpz_t(), upperhalf_cf_terms.continuants.q_k1.get_mpz_t());

  if (upper_half_size > Params::ThresholdAddSubInParallel)	// Subtract in parallel
  {
    std::thread th(mpz_sub, correction_num.get_mpz_t(), c_num1.get_mpz_t(), c_num2.get_mpz_t());
    mpz_sub(correction_den.get_mpz_t(), c_den1.get_mpz_t(), c_den2.get_mpz_t());
    th.join();
  }
  else
  {
    mpz_sub(correction_num.get_mpz_t(), c_num1.get_mpz_t(), c_num2.get_mpz_t());
    mpz_sub(correction_den.get_mpz_t(), c_den1.get_mpz_t(), c_den2.get_mpz_t());
  }

	if (mpz_sgn(correction_num.get_mpz_t()) < 0 && mpz_sgn(correction_den.get_mpz_t()) < 0)
	{
		mpz_neg(correction_num.get_mpz_t(), correction_num.get_mpz_t());
		mpz_neg(correction_den.get_mpz_t(), correction_den.get_mpz_t());
	}

	// Correct CF of Upper half if correction term is negative
	int corrections = 0;
  if (mpz_sgn(correction_num.get_mpz_t()) < 0 || mpz_sgn(correction_den.get_mpz_t()) < 0 || (correction_num < correction_den))
  {
    corrections = perform_corrections(correction_num, correction_den, upperhalf_cf_terms);

    if (upperhalf_cf_terms.list.empty())
		{
      upperhalf_cf_terms.terms_on_disk -= corrections;
      std::filesystem::resize_file(upperhalf_cf_terms.get_filename(), sizeof(mp_limb_t) * upperhalf_cf_terms.terms_on_disk);
    }

    upperhalf_cf_terms.updated = false;

    assert(correction_num >= 0 && correction_den > 0);

    if (correction_num == 0)
      return upperhalf_cf_terms;
  }

	assert(correction_num >= correction_den);

	//if (!half) {
	//	std::cout << "#" << std::flush;
	//}

	CFTermList lowerhalf_cf_terms;

	if (half)
	{
		shift = std::min(mpz_size(correction_num.get_mpz_t()), mpz_size(correction_den.get_mpz_t()));

		if (     mpz_size(correction_num.get_mpz_t()) < Params::ThresholdLowerHalfReductionSize1)
			shift /= Params::LowerHalfReductionFactor1;
		else if (mpz_size(correction_num.get_mpz_t()) < Params::ThresholdLowerHalfReductionSize2)
			shift /= Params::LowerHalfReductionFactor2;
		else if (mpz_size(correction_num.get_mpz_t()) < Params::ThresholdLowerHalfReductionSize3)
			shift /= Params::LowerHalfReductionFactor3;
		else if (mpz_size(correction_num.get_mpz_t()) < Params::ThresholdLowerHalfReductionSize4)
			shift /= Params::LowerHalfReductionFactor4;
		else shift /= 3;

		shift *= sizeof(mp_limb_t) * 8;

		lowerhalf_cf_terms = reg_cf_expansion(correction_num >> shift, correction_den >> shift, calc_convergents, true);
	}
	else lowerhalf_cf_terms = reg_cf_expansion(correction_num, correction_den, calc_convergents);

	if (lowerhalf_cf_terms.list.size() > Params::CFTermsUseDisk)
		lowerhalf_cf_terms.offload();
	
	if (calc_convergents && !upperhalf_cf_terms.updated)
		upperhalf_cf_terms.update_continuants(corrections == 1);

	upperhalf_cf_terms.append(lowerhalf_cf_terms);
	
	if (!calc_convergents) { upperhalf_cf_terms.updated = false; return upperhalf_cf_terms; }

	if (!upperhalf_cf_terms.updated && !lowerhalf_cf_terms.updated)
	{
		upperhalf_cf_terms.update_continuants();
		return upperhalf_cf_terms;
	}
	else
	{
		//mpz_class p_k, q_k, p_k1, q_k1;

		// TODO Move to Continuant namespace via a combine_continuants() func
  //  multiplication::mul( p_k.get_mpz_t(), upperhalf_cf_terms.p_k.get_mpz_t(), lowerhalf_cf_terms. p_k.get_mpz_t());
  //  multiplication::mul( q_k.get_mpz_t(), upperhalf_cf_terms.q_k.get_mpz_t(), lowerhalf_cf_terms. p_k.get_mpz_t());
  //  multiplication::mul(p_k1.get_mpz_t(), upperhalf_cf_terms.p_k.get_mpz_t(), lowerhalf_cf_terms.p_k1.get_mpz_t());
  //  multiplication::mul(q_k1.get_mpz_t(), upperhalf_cf_terms.q_k.get_mpz_t(), lowerhalf_cf_terms.p_k1.get_mpz_t());

  //  multiplication::addmul( p_k.get_mpz_t(), upperhalf_cf_terms.p_k1.get_mpz_t(), lowerhalf_cf_terms. q_k.get_mpz_t());
  //  multiplication::addmul( q_k.get_mpz_t(), upperhalf_cf_terms.q_k1.get_mpz_t(), lowerhalf_cf_terms. q_k.get_mpz_t());
  //  multiplication::addmul(p_k1.get_mpz_t(), upperhalf_cf_terms.p_k1.get_mpz_t(), lowerhalf_cf_terms.q_k1.get_mpz_t());
  //  multiplication::addmul(q_k1.get_mpz_t(), upperhalf_cf_terms.q_k1.get_mpz_t(), lowerhalf_cf_terms.q_k1.get_mpz_t());

		//upperhalf_cf_terms.p_k  =  p_k;
		//upperhalf_cf_terms.q_k  =  q_k;
		//upperhalf_cf_terms.p_k1 = p_k1;
		//upperhalf_cf_terms.q_k1 = q_k1;
		combine_continuants(upperhalf_cf_terms.continuants, lowerhalf_cf_terms.continuants);
		upperhalf_cf_terms.updated = true;
	}

	return upperhalf_cf_terms;
}

void crunch_reg_cf_expansion(std::string file, size_t terms)
{
	mpq_class frac;
	int digits_base = 0;
	double split_factor = 0;

	while (digits_base < 1 || digits_base > 2) {
		std::cout << "The digits are stored in:\n1 - Decimal\n2 - Hexadecimal\nEnter: ";
		std::cin >> digits_base;
	}

	std::cout << "NOTE: Split factor must be between (0, 1). A higher split factor corresponds\n"
		        << "to lesser RAM usage, but slightly slower performance."
		        << std::endl;
	do {
		std::cout << "Enter the split factor (default = 0.5): ";
		std::cin >> split_factor;
	} while (split_factor <= 0 || split_factor >= 1);

	std::cout << "Initializing multiplication tables...";
	multiplication::init();
	std::cout << " done" << std::endl;

	std::cout << "Reading from file '" << file << "'...";
	read_from_file(&frac, file, digits_base == 1);
	std::cout << " done [frac = " << frac.get_d() << "]" << std::endl;

	double session_time = wall_clock();

	int iteration = 0;
	size_t computed_terms = 0;

	while (frac.get_num() > 0 && frac.get_den() > 0)
	{
		iteration++;
		std::cerr << "Iteration " << iteration << "...";

		double iter_start_time = wall_clock();

		size_t redc = std::min(mpz_size(frac.get_num_mpz_t()), mpz_size(frac.get_den_mpz_t()));
		redc = size_t(redc * split_factor);

		size_t shift = redc * sizeof(mp_limb_t) * 8;

		CFTermList convergents(reg_cf_expansion(frac.get_num() >> shift, frac.get_den() >> shift, true, true));

		if (convergents.list.empty())
			convergents.onload();

		if (!convergents.updated) convergents.update_continuants();

		// Determine correction fraction
		{
			// Numerator = b * p_{k-1} - a * q_{k-1}
			// Denomenator = a * q_k - b * p_k
			mpz_class c_num, c_den;

			multiplication::   mul(c_num.get_mpz_t(), frac.get_den_mpz_t(), convergents.continuants.p_k1.get_mpz_t());
			multiplication::submul(c_num.get_mpz_t(), frac.get_num_mpz_t(), convergents.continuants.q_k1.get_mpz_t());

			multiplication::   mul(c_den.get_mpz_t(), frac.get_num_mpz_t(), convergents.continuants.q_k.get_mpz_t());
			multiplication::submul(c_den.get_mpz_t(), frac.get_den_mpz_t(), convergents.continuants.p_k.get_mpz_t());

			if (mpz_sgn(c_num.get_mpz_t()) < 0 && mpz_sgn(c_den.get_mpz_t()) < 0)
			{
				mpz_neg(c_num.get_mpz_t(), c_num.get_mpz_t());
				mpz_neg(c_den.get_mpz_t(), c_den.get_mpz_t());
			}

			frac.get_num() = c_num;
			frac.get_den() = c_den;
		}

		// Correct CF of Upper half if correction term is negative
		size_t corrections = 0;
		if (sgn(frac.get_num()) < 0 || sgn(frac.get_den()) < 0 || (frac.get_num() < frac.get_den()))
		{
			corrections = perform_corrections(frac.get_num(), frac.get_den(), convergents/*.list*/);

			convergents.updated = false;

			assert(frac.get_num() > 0);
			assert(frac.get_den() > 0);
		}

		assert(frac.get_num() >= frac.get_den());

		computed_terms += convergents.list.size();

		double iter_end_time = wall_clock() - iter_start_time;

		std::cerr.flags(std::cerr.fixed);
		std::cerr.precision(2);
		std::cerr << "\t" << corrections << " corrections\t"
			<< iter_end_time << "ms\t" << computed_terms << "/" << terms
			<< "\t(" << (100. * computed_terms) / double(terms) << "%)\t"
			<< int(convergents.list.size() / iter_end_time) << "K Terms/s\n";

		output_cf_terms_list(convergents.list, std::string("iteration") + std::to_string(iteration) + std::string(".txt"), false);

		if (computed_terms >= terms) break;
	}

	session_time = wall_clock() - session_time;

	std::cout << "Computed " << computed_terms << " terms in " << iteration
		<< " iterations and " << session_time << "s" << std::endl;

	std::cin.get();
}

void crunch_reg_cf_expansion_on_disk(std::string file, size_t terms, size_t nthreads)
{
	bool resuming_computation = file == "";

	if (!resuming_computation)
		std::cerr << "Computing " << terms << " terms from file '" << file << "' using "
		<< (double)disk_mpz::SplitSize / (1024 * 1024) << "MB/file" << std::endl;
	else
		std::cerr << "Resuming computation of " << terms << " terms using " << (double)disk_mpz::SplitSize / (1024. * 1024.) <<
		"MB/file" << std::endl;

	int iteration = 1;
	size_t computed_terms = 0;

	if (resuming_computation)
	{
		// Check iteration number and computed terms
		while (true)
		{
			std::string iter_f("iterations\\iteration" + std::to_string(iteration) + ".txt");

			if (file_exists(iter_f))
			{
				iteration++;
				computed_terms += count_lines(iter_f);
			}
			else break;
		}
	}

	double session_time = wall_clock();

	disk_mpq frac("frac", file);

	std::cout << std::endl;

	std::cout << "Initializing multiplication tables...";
	multiplication::init();
	std::cout << " done" << std::endl;

	while (true)
	{
		std::cerr << "Iteration " << iteration << "...\n---\nComputing terms...\t";

		double iter_start_time = wall_clock();

		// Get the two most significant limbs of the fraction as numerator and denominator
		mpz_class num = frac.get_num().get_top2_mpz();
		mpz_class den = frac.get_den().get_top2_mpz();

		// Shift the numerator and denominator so as to convert them into standard size
		size_t shift = mpz_sizeinbase(num.get_mpz_t(), 2) - (disk_mpz::SplitSize * 8);
		num >>= shift;
		den >>= shift;

		// Calculate regular continued fraction of the top two terms of the fraction
		CFTermList cf_expansion (reg_cf_expansion(num, den, true, true));

		cf_expansion.onload();	// Load into RAM

		for (int i = 0; i < Params::TruncateConvergents; ++i)
			cf_expansion.list.pop_back();
		cf_expansion.update_continuants();

		std::cerr << cf_expansion.list.size() << " terms | " << int(wall_clock() - iter_start_time) << "ms | ";
		std::cerr << cf_expansion.list.size() / int(wall_clock() - iter_start_time) << "KT/s | ";
		std::cerr << "largest: " << *std::max_element(cf_expansion.list.begin(), cf_expansion.list.end()) << std::endl;
		std::cerr << std::endl;
		
		{
			double t, end_t;
			std::cerr.precision(4);
			std::cerr << "Calculating correction fraction" << std::endl;

			// Determine correction fraction
			size_t nfiles = frac.get_num().files() + frac.get_den().files();

			std::cerr << "numerator...\t";

			t = wall_clock();

			// correction numerator = den * p_k1 - num * q_k1
			disk_mpz c_num(disk_mpz::cross_mul_sub("corr_num", frac.get_den(), cf_expansion.continuants.p_k1, frac.get_num(), cf_expansion.continuants.q_k1, nthreads));
			end_t = wall_clock() - t;
			
			std::cerr<< int(end_t) << "ms (" << (disk_mpz::SplitSize * nfiles) / (1024. * 1.024 * end_t) << "MB/s)" << std::endl;

			std::cerr << "denominator...\t";

			t = wall_clock();
			// correction denominator = num * q_k - den * p_k
			disk_mpz c_den(disk_mpz::cross_mul_sub("corr_den", frac.get_num(), cf_expansion.continuants.q_k, frac.get_den(), cf_expansion.continuants.p_k, nthreads));
			end_t = wall_clock() - t;

			std::cerr << int(end_t) << "ms (" << (disk_mpz::SplitSize * nfiles) / (1024. * 1.024 * end_t) << "MB/s)" << std::endl;

			if (c_num.sign() < 0 && c_den.sign() < 0)
			{
				c_num.set_pos();
				c_den.set_pos();
			}
			else if (c_num.sign() < 0 || c_den.sign() < 0)
				throw "Correction fraction negative";

			std::cerr << "moving...\t";

			t = wall_clock();
			disk_mpz_move(frac.get_num(), c_num);
			disk_mpz_move(frac.get_den(), c_den);
			end_t = wall_clock();

			std::cerr << int(end_t - t) << "ms" << std::endl;
			std::cerr << std::endl;
		}

		// No corrections to perform as correction fraction is almost certainly positive

		computed_terms += cf_expansion.list.size();
		double iter_end_time = wall_clock() - iter_start_time;

		output_cf_terms_list(cf_expansion.list, std::string("iteration") + std::to_string(iteration) + std::string(".txt"), false);

		std::cerr.flags(std::cerr.fixed);
		std::cerr.precision(2);

//  std::cerr << "Scratch memory usage for multiplication: " << multiplication::scratch_mem_size / (1024.*1024.) << "MB" << std::endl;
		std::cerr << "Completed " << (100. * computed_terms) / double(terms) << "% in " << (int)iter_end_time << "ms ";
		std::cerr << "(" << int(cf_expansion.list.size() / iter_end_time) << "K Terms/s)" << std::endl;
		std::cerr << std::endl;

		if (computed_terms >= terms) break;

		iteration++;
	}

	session_time = wall_clock() - session_time;

	std::cout << "Computed " << computed_terms << " terms in " << iteration << " iterations and " << session_time << "ms ";
	std::cout << "(" << computed_terms / session_time << "K Terms/s)" << std::endl;

	std::cin.get();
}

void CFTermList::update_continuants(bool single_term)
{
	bool keep_loaded = false;
	if (list.size() == 0)
		onload();
	else keep_loaded = true;

	if (list.size() > 2 && list.size() < Params::ThresholdUseBasicContProc) {
		auto lo = continuant_pair_right(0, list.size() - 1, list);
		auto hi = continuant_pair_right(1, list.size() - 1, list);
		continuants.p_k = lo.first;
		continuants.p_k1 = lo.second;
		continuants.q_k = hi.first;
		continuants.q_k1 = hi.second;
	}
	else {
		auto mid = list.size() / 2;
		ContinuantCache::clear();

		if (list.size() > 1)
		{
			if (!single_term) ContinuantCache::build_cache_indices(0, list.size() - 1, mid);
			ContinuantCache::build_cache_indices(0, list.size() - 2, mid);
			if (!single_term) ContinuantCache::build_cache_indices(1, list.size() - 1, mid);
			ContinuantCache::build_cache_indices(1, list.size() - 2, mid);
		}

		if (single_term) continuants.p_k = continuants.p_k1;
		else             continuants.p_k = list.size() == 1 ? 1 : continuant(0, list.size() - 1, list, mid);

		continuants.p_k1 = list.size() == 1 ? 1 : continuant(0, list.size() - 2, list, mid);

		if (single_term) continuants.q_k = continuants.q_k1;
		else             continuants.q_k = list.size() == 1 ? 1 : continuant(1, list.size() - 1, list, mid);

		continuants.q_k1 = list.size() == 1 ? 0 : list.size() == 2 ? 1 : continuant(1, list.size() - 2, list, mid);

		ContinuantCache::clear();
	}

	updated = true;

	if (list.size() > Params::CFTermsUseDisk && !keep_loaded)
		offload();
}

// offload() moves a list from RAM to disk
void CFTermList::offload()
{
	static size_t cf_terms_file_idx = 1;
	terms_on_disk = list.size();

	if (idx == 0)
		idx = cf_terms_file_idx++;

	std::ofstream cfterms_file(get_filename(), std::ios::binary | std::ios::out);

	cfterms_file.write((const char*) list.data(), sizeof(mpir_ui) * list.size());

	cfterms_file.close();
	list.clear();
}

// Load data from disk to RAM
void CFTermList::onload()
{
	list.resize(terms_on_disk);

	std::ifstream cfterms_file(get_filename(), std::ios::binary | std::ios::in);
	cfterms_file.read((char*)list.data(), terms_on_disk * sizeof(mp_limb_t));

	cfterms_file.close();

	std::filesystem::remove(get_filename());
}

// Append more terms from another list
// If either of the lists is on disk, then appending
// takes place on disk. Otherwise the list is kept on RAM
void CFTermList::append(const CFTermList& to_append)
{
	auto append_to_file = [&]()
	{
		std::filesystem::resize_file(get_filename(), sizeof(mp_limb_t) * terms_on_disk);
		std::ofstream outFile(get_filename(), std::ios::app | std::ios::binary);

		if (to_append.terms_on_disk > 0)
		{
			std::ifstream inFile(to_append.get_filename(), std::ios::binary | std::ios::in);

			outFile << inFile.rdbuf();
			inFile.close();

			terms_on_disk += to_append.terms_on_disk;

			std::filesystem::remove(to_append.get_filename());
		}
		else
		{
			outFile.write((const char*)to_append.list.data(), sizeof(mp_limb_t) * to_append.list.size());
			terms_on_disk += to_append.list.size();
		}
		outFile.close();
	};

	if (terms_on_disk > 0) // If already offloaded
		append_to_file();
	else
	{
		if (to_append.terms_on_disk == 0)
			list.insert(list.end(), to_append.list.begin(), to_append.list.end());
		else // Offload if the other list is already offloaded
		{
			offload();
			append_to_file();
		}
	}
}