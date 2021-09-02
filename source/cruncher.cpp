#include "cruncher.h"
#include "continuant.h"
#include "io.h"
#include "utils.h"
#include "disk_mpz.h"
#include <cassert>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <thread>
#include <filesystem>

void basecase_reg_cf_terms(const mpz_class& num, const mpz_class& den, CFTerms& out, bool half) {
	mpz_class r, N(num), D(den), s;
	mpz_class* n = &N, * d = &D;

	if (half) s = sqrt(num);

	while (*d != 0) {
		mpz_tdiv_q(r.get_mpz_t(), n->get_mpz_t(), d->get_mpz_t());

		if (!r.fits_ui_p()) throw "Unusually large convergent";

		out.push_back(r.get_ui());
		mpz_submul(n->get_mpz_t(), d->get_mpz_t(), r.get_mpz_t());

		if (half && *n < s) break;

		std::swap(n, d);
	}
}

int perform_correction(mpz_class& num, mpz_class& den, CFTermList& c) {
	int corrections = 0;
	while (mpz_sgn(num.get_mpz_t()) < 0 || mpz_sgn(den.get_mpz_t()) < 0 || (num < den)) {
		corrections++;
		mpz_swap(num.get_mpz_t(), den.get_mpz_t());

		if (mpz_sgn(den.get_mpz_t()) < 0) {
			assert(num > 0);
			mpz_neg(num.get_mpz_t(), num.get_mpz_t());
			mpz_neg(den.get_mpz_t(), den.get_mpz_t());
		}

		if (c.list.empty()) {
			uint64_t lastterm;
			std::ifstream cfterms_file(c.get_filename(), std::ios::in | std::ios::binary);
			cfterms_file.seekg((c.terms_on_disk - corrections) * sizeof(uint64_t));
			cfterms_file.read((char*)&lastterm, sizeof(uint64_t));
			mpz_addmul_ui(num.get_mpz_t(), den.get_mpz_t(), lastterm);
		}
		else {
			mpz_addmul_ui(num.get_mpz_t(), den.get_mpz_t(), c.list.back());
			c.list.pop_back();

			if (c.list.empty())
				return corrections;
		}
	}
	return corrections;
}

CFTermList reg_cf_terms(const mpz_class& num, const mpz_class& den, bool calc_convergents, bool half) {
	assert(num >= 0);
	assert(den > 0);

	size_t reduction_size = std::min(mpz_size(num.get_mpz_t()), mpz_size(den.get_mpz_t())) / 2;

	if (reduction_size < Params::ThresholdCFDivide) {
		CFTermList out;
		out.list.reserve(Params::DefaultConvAllocSize);
		basecase_reg_cf_terms(num, den, out.list, half);
		if (calc_convergents)
			out.update();
		return out;
	}

	// Take CF of upper half
	CFTermList upperhalf_cf_terms(reg_cf_terms(num >> (reduction_size * sizeof(mp_limb_t) * 8),
		den >> (reduction_size * sizeof(mp_limb_t) * 8), true, true));

	assert(upperhalf_cf_terms.terms_on_disk > 0 || !upperhalf_cf_terms.list.empty());
	
	if (upperhalf_cf_terms.list.size() > Params::CFTermsUseDisk)
		upperhalf_cf_terms.offload();
	
	// Determine correction fraction
	// Numerator = b * p_{k-1} - a * q_{k-1}
	// Denomenator = a * q_k - b * p_k
	assert(upperhalf_cf_terms.updated);

	// Launch a new thread to calculate denom
	mpz_class correction_num(0), correction_den(0);

	if (reduction_size * 2 > Params::CFTermsUseParallelMultThreshold1) {
		// correction_den = c_den1 - c_den2 (same for correction_num)
		mpz_class c_den1, c_den2, c_num1, c_num2;

		std::thread th1(mpz_mul, c_den1.get_mpz_t(), num.get_mpz_t(), upperhalf_cf_terms.q_k.get_mpz_t());
		std::thread th2(mpz_mul, c_den2.get_mpz_t(), den.get_mpz_t(), upperhalf_cf_terms.p_k.get_mpz_t());
		std::thread th3(mpz_mul, c_num1.get_mpz_t(), den.get_mpz_t(), upperhalf_cf_terms.p_k1.get_mpz_t());
		mpz_mul(c_num2.get_mpz_t(), num.get_mpz_t(), upperhalf_cf_terms.q_k1.get_mpz_t());

		th3.join();
		mpz_sub(correction_num.get_mpz_t(), c_num1.get_mpz_t(), c_num2.get_mpz_t());

		th1.join();
		th2.join();
		mpz_sub(correction_den.get_mpz_t(), c_den1.get_mpz_t(), c_den2.get_mpz_t());
	}
	else {
		mpz_addmul(correction_den.get_mpz_t(), num.get_mpz_t(), upperhalf_cf_terms.q_k.get_mpz_t());
		mpz_submul(correction_den.get_mpz_t(), den.get_mpz_t(), upperhalf_cf_terms.p_k.get_mpz_t());

		mpz_addmul(correction_num.get_mpz_t(), den.get_mpz_t(), upperhalf_cf_terms.p_k1.get_mpz_t());
		mpz_submul(correction_num.get_mpz_t(), num.get_mpz_t(), upperhalf_cf_terms.q_k1.get_mpz_t());
	}

	if (mpz_sgn(correction_num.get_mpz_t()) < 0 && mpz_sgn(correction_den.get_mpz_t()) < 0) {
		mpz_neg(correction_num.get_mpz_t(), correction_num.get_mpz_t());
		mpz_neg(correction_den.get_mpz_t(), correction_den.get_mpz_t());
	}

	// Correct CF of Upper half if correction term is negative
	int corrections = 0;
	if (mpz_sgn(correction_num.get_mpz_t()) < 0
		|| mpz_sgn(correction_den.get_mpz_t()) < 0
		|| (correction_num < correction_den)) {
		corrections = perform_correction(correction_num, correction_den, upperhalf_cf_terms);

		
		if (upperhalf_cf_terms.list.empty()) {
			upperhalf_cf_terms.terms_on_disk -= corrections;
			std::filesystem::resize_file(upperhalf_cf_terms.get_filename(),
				sizeof(uint64_t) * upperhalf_cf_terms.terms_on_disk);
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

	if (half) {
		reduction_size = std::min(mpz_size(correction_num.get_mpz_t()),
			mpz_size(correction_den.get_mpz_t()));	// Reduction terms_on_disk

		if (mpz_size(correction_num.get_mpz_t()) < Params::ThresholdReduc1)
			reduction_size /= 2;
		else if (mpz_size(correction_num.get_mpz_t()) < Params::ThresholdReduc2)
			reduction_size /= 2.8;
		else if (mpz_size(correction_num.get_mpz_t()) < Params::ThresholdReduc3)
			reduction_size /= 3;
		else if (mpz_size(correction_num.get_mpz_t()) < Params::ThresholdReduc4)
			reduction_size /= 3;
		else reduction_size /= 3;

		lowerhalf_cf_terms = reg_cf_terms(correction_num >> (reduction_size * sizeof(mp_limb_t) * 8),
			correction_den >> (reduction_size * sizeof(mp_limb_t) * 8), calc_convergents, true);
	}
	else lowerhalf_cf_terms = reg_cf_terms(correction_num, correction_den, calc_convergents);

	
	if (lowerhalf_cf_terms.list.size() > Params::CFTermsUseDisk)
		lowerhalf_cf_terms.offload();
	

	if (calc_convergents && !upperhalf_cf_terms.updated)
		upperhalf_cf_terms.update(corrections == 1);

	upperhalf_cf_terms.append(lowerhalf_cf_terms);
	

	if (!calc_convergents) { upperhalf_cf_terms.updated = false; return upperhalf_cf_terms; }

	if (!upperhalf_cf_terms.updated && !lowerhalf_cf_terms.updated) {
		upperhalf_cf_terms.update();
		return upperhalf_cf_terms;
	}
	else {
		mpz_class p_k, q_k, p_k1, q_k1;

		if (reduction_size * 2 <= Params::CFTermsUseParallelMultThreshold2) {
			p_k = upperhalf_cf_terms.p_k * lowerhalf_cf_terms.p_k
				+ upperhalf_cf_terms.p_k1 * lowerhalf_cf_terms.q_k;
			q_k = upperhalf_cf_terms.q_k * lowerhalf_cf_terms.p_k
				+ upperhalf_cf_terms.q_k1 * lowerhalf_cf_terms.q_k;
			p_k1 = upperhalf_cf_terms.p_k * lowerhalf_cf_terms.p_k1
				+ upperhalf_cf_terms.p_k1 * lowerhalf_cf_terms.q_k1;
			q_k1 = upperhalf_cf_terms.q_k * lowerhalf_cf_terms.p_k1
				+ upperhalf_cf_terms.q_k1 * lowerhalf_cf_terms.q_k1;
		}
		else {
			std::thread th1(mpz_mul, p_k.get_mpz_t(),
				upperhalf_cf_terms.p_k.get_mpz_t(), lowerhalf_cf_terms.p_k.get_mpz_t());
			std::thread th2(mpz_mul, q_k.get_mpz_t(),
				upperhalf_cf_terms.q_k.get_mpz_t(), lowerhalf_cf_terms.p_k.get_mpz_t());
			std::thread th3(mpz_mul, p_k1.get_mpz_t(),
				upperhalf_cf_terms.p_k.get_mpz_t(), lowerhalf_cf_terms.p_k1.get_mpz_t());
			mpz_mul(q_k1.get_mpz_t(),
				upperhalf_cf_terms.q_k.get_mpz_t(), lowerhalf_cf_terms.p_k1.get_mpz_t());

			th1.join();
			th1 = std::thread(mpz_addmul, p_k.get_mpz_t(),
				upperhalf_cf_terms.p_k1.get_mpz_t(), lowerhalf_cf_terms.q_k.get_mpz_t());
			th2.join();
			th2 = std::thread(mpz_addmul, q_k.get_mpz_t(),
				upperhalf_cf_terms.q_k1.get_mpz_t(), lowerhalf_cf_terms.q_k.get_mpz_t());
			th3.join();
			th3 = std::thread(mpz_addmul, p_k1.get_mpz_t(),
				upperhalf_cf_terms.p_k1.get_mpz_t(), lowerhalf_cf_terms.q_k1.get_mpz_t());
			mpz_addmul(q_k1.get_mpz_t(),
				upperhalf_cf_terms.q_k1.get_mpz_t(), lowerhalf_cf_terms.q_k1.get_mpz_t());

			th1.join();
			th2.join();
			th3.join();
		}

		upperhalf_cf_terms.p_k = p_k;
		upperhalf_cf_terms.q_k = q_k;
		upperhalf_cf_terms.p_k1 = p_k1;
		upperhalf_cf_terms.q_k1 = q_k1;
		upperhalf_cf_terms.updated = true;
	}
	return upperhalf_cf_terms;
}

void crunch_reg_cf_terms(std::string file, size_t terms) {
	int c = 0;
	while (c < 1 || c > 2) {
		std::cout << "The digits are stored in:\n1 - Decimal\n2 - Hexadecimal\nEnter: ";
		std::cin >> c;
	}

	int target_iterations = -1;
	do {
		std::cout << "Larger # of target iterations = slow speed, but less RAM consumed.\n";
		std::cout << "Enter the number of target iterations (very imprecise, enter 0 to use default): ";
		std::cin >> target_iterations;

		if (target_iterations == 0) target_iterations = 10;
	} while (target_iterations == -1);

	mpq_class frac;
	read_from_file(&frac, file, c == 1);
	int iter = 0;
	size_t computed_terms = 0;
	double session_time = wall_clock();

	do {
		iter++;
		std::cerr << "Iteration " << iter << "...";

		double iter_start_time = wall_clock();

		size_t redc = std::min(mpz_size(frac.get_num_mpz_t()), mpz_size(frac.get_den_mpz_t()));
		redc = size_t(redc / (1. + 1. / target_iterations));	// Reduction terms_on_disk

		CFTermList convergents(reg_cf_terms(frac.get_num() >> (redc * sizeof(mp_limb_t) * 8),
			frac.get_den() >> (redc * sizeof(mp_limb_t) * 8), true, true));

		if (convergents.list.empty())
			convergents.onload();

		if (!convergents.updated) convergents.update();

		{
			// Determine correction fraction
			// Numerator = b * p_{k-1} - a * q_{k-1}
			// Denomenator = a * q_k - b * p_k
			mpz_class c_num;
			mpz_realloc(c_num.get_mpz_t(), mpz_size(frac.get_den_mpz_t()) + mpz_size(convergents.p_k1.get_mpz_t()));
			c_num = (frac.get_den() * convergents.p_k1);
			mpz_submul(c_num.get_mpz_t(), frac.get_num_mpz_t(), convergents.q_k1.get_mpz_t());
			mpz_class c_den(frac.get_num() * convergents.q_k);
			mpz_submul(c_den.get_mpz_t(), frac.get_den_mpz_t(), convergents.p_k.get_mpz_t());

			if (mpz_sgn(c_num.get_mpz_t()) < 0 && mpz_sgn(c_den.get_mpz_t()) < 0) {
				mpz_neg(c_num.get_mpz_t(), c_num.get_mpz_t());
				mpz_neg(c_den.get_mpz_t(), c_den.get_mpz_t());
			}

			frac.get_num() = c_num;
			frac.get_den() = c_den;
		}

		// Correct CF of Upper half if correction term is negative
		size_t corrections = 0;
		if (sgn(frac.get_num()) < 0 || sgn(frac.get_den()) < 0 || (frac.get_num() < frac.get_den())) {
			corrections = perform_correction(frac.get_num(), frac.get_den(), convergents/*.list*/);

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

		output_cf_terms_list(convergents.list, std::string("iteration") + to_str(iter) + std::string(".txt"), false);

		if (computed_terms >= terms)
			break;
	} while (frac.get_num() > 0 && frac.get_den() > 0);

	session_time = wall_clock() - session_time;

	std::cout << "Computed " << computed_terms << " terms in " << iter
		<< " iterations and " << session_time << "s" << std::endl;

	std::cin.get();
}

void crunch_reg_cf_terms_on_disk(std::string file, size_t terms, size_t bytes_per_file, size_t nthreads) {
	if (file != "")
		std::cerr << "Computing " << terms << " terms from file '" << file << "' using "
		<< (double)bytes_per_file / (1024 * 1024) << "MB/file" << std::endl;
	else
		std::cerr << "Resuming computation of " << terms << " terms using " << (double)bytes_per_file / (1024. * 1024.) <<
		"MB/file" << std::endl;

	int iter = 1;
	size_t computed_terms = 0;

	if (file == "") {	// Check iteration number and computed terms
		while (true) {
			std::string iter_f("iterations\\iteration" + std::to_string(iter) + ".txt");

			if (file_exists(iter_f)) {
				iter++;
				computed_terms += count_lines(iter_f);
			}
			else break;
		}
	}
	double session_time = wall_clock();

	disk_mpq frac("frac", file, bytes_per_file);

	do {
		std::cerr << "Iteration " << iter << "...\nComputing terms...";

		double iter_start_time = wall_clock();

		CFTermList convergents
		(reg_cf_terms(frac.get_num().get_top2_mpz(), frac.get_den().get_top2_mpz(), true, true));

		convergents.onload();

		for (int i = 0; i < Params::TruncateConvergents; ++i)
			convergents.list.pop_back();
		convergents.update();

		std::cerr << "\t" << int(wall_clock() - iter_start_time) << "ms for\t" << convergents.list.size() << " terms \t("
			<< convergents.list.size() / int(wall_clock() - iter_start_time) << "KT/s)    largest: "
			<< *std::max_element(convergents.list.begin(), convergents.list.end()) << std::endl;

		{
			std::cerr << "Multiplying...";

			// Determine correction fraction
			size_t nfiles = frac.get_num().files() + frac.get_den().files();

			double t = wall_clock();
			disk_mpz c_num(disk_mpz::cross_mult_sub("corr_num",
				frac.get_den(), convergents.p_k1, frac.get_num(), convergents.q_k1,
					bytes_per_file, nthreads));

			std::cerr.precision(4);
			double end_t = wall_clock() - t;
			std::cerr << "\t" << int(end_t) << "ms (den, \t"
				<< (bytes_per_file * nfiles) / (1024. * 1.024 * end_t) << "MB/s)...";

			t = wall_clock();
			disk_mpz c_den(disk_mpz::cross_mult_sub("corr_den",
				frac.get_num(), convergents.q_k, frac.get_den(), convergents.p_k,
				bytes_per_file, nthreads));
			end_t = wall_clock() - t;

			std::cerr << "\t" << int(end_t) << "ms (num, \t"
				<< (bytes_per_file * nfiles) / (1024. * 1.024 * end_t) << "MB/s)...";

			if (c_num.sign() < 0 && c_den.sign() < 0) {
				c_num.set_pos(); c_den.set_pos();
			}
			else if (c_num.sign() < 0 || c_den.sign() < 0) {
				throw "Correction fraction negative";
			}

			t = wall_clock();
			disk_mpz_move(frac.get_num(), c_num);
			disk_mpz_move(frac.get_den(), c_den);
			std::cerr << "\t" << int(wall_clock() - t) << "ms (move)..." << std::endl;;
		}

		// No corrections to perform as correction fraction is almost certainly positive

		computed_terms += convergents.list.size();
		double iter_end_time = wall_clock() - iter_start_time;

		output_cf_terms_list(convergents.list, std::string("iteration") + to_str(iter) + std::string(".txt"), false);

		std::cerr.flags(std::cerr.fixed);
		std::cerr.precision(2);
		std::cerr << "Completed \t" << (100. * computed_terms) / double(terms) << "% in \t"
			<< (int)iter_end_time << "ms\t(" 
			<< int(convergents.list.size() / iter_end_time) << "K Terms/s)\n\n";

		if (computed_terms >= terms)
			break;
		iter++;
	} while (true);

	session_time = wall_clock() - session_time;

	std::cout << "Computed " << computed_terms << " terms in " << iter
		<< " iterations and " << session_time << "s" << std::endl;

	std::cin.get();
}


void CFTermList::update(bool single_term) {
	bool keep_loaded = false;
	if (list.size() == 0)
		onload();
	else keep_loaded = true;

	if (list.size() > 2 && list.size() < Params::ThresholdUseBasicContProc) {
		auto lo = continuant_pair_right(0, list.size() - 1, list);
		auto hi = continuant_pair_right(1, list.size() - 1, list);
		p_k = lo.first;
		p_k1 = lo.second;
		q_k = hi.first;
		q_k1 = hi.second;
	}
	else {
		auto mid = list.size() / 2;
		ContinuantCache::clear();

		//--
		if (!single_term) ContinuantCache::dummy_run(0, list.size() - 1, mid);
		ContinuantCache::dummy_run(0, list.size() - 2, mid);
		if (!single_term) ContinuantCache::dummy_run(1, list.size() - 1, mid);
		ContinuantCache::dummy_run(1, list.size() - 2, mid);
		//--

		if (single_term) p_k = p_k1;
		else p_k = list.size() == 1 ? 1 : continuant(0, list.size() - 1, list, mid);
		p_k1 = list.size() == 1 ? 1 : continuant(0, list.size() - 2, list, mid);

		if (single_term) q_k = q_k1;
		else q_k = list.size() == 1 ? 1 : continuant(1, list.size() - 1, list, mid);
		q_k1 = list.size() == 1 ? 0 : list.size() == 2 ? 1 : continuant(1, list.size() - 2, list, mid);

		ContinuantCache::clear();
	}

	updated = true;

	if (list.size() > Params::CFTermsUseDisk && !keep_loaded)
		offload();
}

void CFTermList::offload()
{
	static size_t cf_terms_file_idx = 1;
	terms_on_disk = list.size();
	if (idx == 0)
		idx = cf_terms_file_idx++;

	std::ofstream cfterms_file(get_filename(),
		std::ios::binary | std::ios::out);

	cfterms_file.write((const char*) list.data(),
		sizeof(mpir_ui) * list.size());

	cfterms_file.close();
	list.clear();
}

void CFTermList::onload() {
	// Load data from disk to RAM
	list.resize(terms_on_disk);

	std::ifstream cfterms_file(get_filename(), std::ios::binary | std::ios::in);
	cfterms_file.read((char*)list.data(), terms_on_disk * sizeof(uint64_t));

	cfterms_file.close();

	std::filesystem::remove(get_filename());
}

void CFTermList::append(const CFTermList& to_append)
{
	auto append_to_file = [&]() {
		std::filesystem::resize_file(get_filename(), sizeof(uint64_t) * terms_on_disk);
		std::ofstream outFile(get_filename(), std::ios::app | std::ios::binary);
		if (to_append.terms_on_disk > 0) {
			std::ifstream inFile(to_append.get_filename(), std::ios::binary | std::ios::in);
			//outFile.seekp();
			outFile << inFile.rdbuf();
			inFile.close();
			terms_on_disk += to_append.terms_on_disk;

			std::filesystem::remove(to_append.get_filename());
		}
		else {
			outFile.write((const char*)to_append.list.data(),
				sizeof(uint64_t) * to_append.list.size());
			terms_on_disk += to_append.list.size();
		}
		outFile.close();
	};

	if (terms_on_disk > 0) {	// If already offloaded
		append_to_file();
	}
	else {
		if (to_append.terms_on_disk == 0)
			list.insert(list.end(), to_append.list.begin(), to_append.list.end());
		else {	// Offload if the other list is already offloaded
			offload();
			append_to_file();
		}
	}
}