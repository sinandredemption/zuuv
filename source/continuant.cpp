#include "continuant.h"
#include <cassert>
#include <thread>
#include <mutex>
#include <atomic>
#include <future>

namespace ContinuantCache {
	std::map<std::pair<size_t, size_t>, mpz_class> table;
	std::map<std::pair<size_t, size_t>, int> table_hit;

	std::mutex mu;
	std::atomic<int> available_threads;
}

void ContinuantCache::dummy_run(size_t s, size_t t, size_t mid) {
	if (table_hit.find({ s,t }) != table_hit.end()) {
		table_hit[{s, t}]++;
		return;
	}

	table_hit[{s, t}] = 0;

	if (t - s <= Params::ThresholdUseBasicContProc)
		return;
	else {
		dummy_run(s, mid, (s + mid) / 2);
		dummy_run(mid + 1, t, (mid + t) / 2);
		dummy_run(s, mid - 1, (s + mid) / 2),
		dummy_run(mid + 2, t, (mid + t) / 2);
	}
}

mpz_class ContinuantCache::single_threaded_continuant(size_t s, size_t t, size_t mid, const CFTerms& terms) {
	if (table.find({ s,t }) != table.end()) {
		table_hit[{s, t}]--;

		if (table_hit[{s, t}] == 0) {
			mpz_class ret(table[{s, t}]);
			table.erase({ s, t });
			return ret;
		}

		return table[{s, t}];
	}

	mpz_class r;
	auto& ret = table_hit[{s, t}] > 0 ? table[{s, t}] : r;

	if (t - s <= Params::ThresholdUseBasicContProc)
		ret = basic_continuant(s, t, terms);
	else {
		ret = single_threaded_continuant(s, mid, (s + mid) / 2, terms)
			* single_threaded_continuant(mid + 1, t, (mid + t) / 2, terms);
		mpz_addmul(ret.get_mpz_t(),
			single_threaded_continuant(s, mid - 1, (s + mid) / 2, terms).get_mpz_t(),
			single_threaded_continuant(mid + 2, t, (mid + t) / 2, terms).get_mpz_t());
	}

	return ret;
}

mpz_class ContinuantCache::multi_threaded_continuant(size_t s, size_t t, size_t mid, const CFTerms& terms) {
	mu.lock();
	if (table.find({ s,t }) != table.end()) {
		mpz_class ret(table[{s, t}]);

		table_hit[{s, t}]--;

		if (table_hit[{s, t}] == 0) {
			table.erase({ s,t });
		}

		mu.unlock();
		return ret;
	}
	mu.unlock();

	mpz_class ret(0);

	if (t - s <= Params::ThresholdUseBasicContProc) {
		ret = basic_continuant(s, t, terms);
	}
	else if (t - s < Params::ContinuantUseParallel || available_threads <= 0) {
		mpz_addmul(ret.get_mpz_t(),
			multi_threaded_continuant(s, mid, (s + mid) / 2, terms).get_mpz_t(),
			multi_threaded_continuant(mid + 1, t, (mid + t) / 2, terms).get_mpz_t());
		mpz_addmul(ret.get_mpz_t(),
			multi_threaded_continuant(s, mid - 1, (s + mid) / 2, terms).get_mpz_t(),
			multi_threaded_continuant(mid + 2, t, (mid + t) / 2, terms).get_mpz_t());
	}
	else {
		mpz_class m[4];
		std::thread th;
		auto new_thread_func = [&]() {
			m[1] = multi_threaded_continuant(mid + 1, t, (mid + t) / 2, terms);
			m[3] = multi_threaded_continuant(mid + 2, t, (mid + t) / 2, terms);
		};

		if (available_threads > 0) {
			available_threads--;
			th = std::thread(new_thread_func);
		}
		else new_thread_func();

		m[0] = multi_threaded_continuant(s, mid, (s + mid) / 2, terms);
		m[2] = multi_threaded_continuant(s, mid - 1, (s + mid) / 2, terms);

		if (th.joinable()) {
			th.join();
			available_threads++;
		}

		mpz_class ret1, ret2;
		if (available_threads) {
			available_threads--;
			th = std::thread(mpz_addmul, ret1.get_mpz_t(), m[0].get_mpz_t(), m[1].get_mpz_t());
		}
		else mpz_addmul(ret1.get_mpz_t(), m[0].get_mpz_t(), m[1].get_mpz_t());
		
		mpz_addmul(ret2.get_mpz_t(), m[2].get_mpz_t(), m[3].get_mpz_t());

		if (th.joinable()) {
			th.join();
			available_threads++;
		}

		mpz_add(ret.get_mpz_t(), ret1.get_mpz_t(), ret2.get_mpz_t());
	}

	mu.lock();
	if (table_hit[{s, t}] > 0) {
		table[{s, t}] = ret;
	}
	mu.unlock();

	return ret;
}

mpz_class continuant(size_t s, size_t t, const CFTerms& terms, size_t split_point) {
	auto mid = split_point == 0 ? (s+t)/2 : split_point;

	using namespace ContinuantCache;

	return multi_threaded_continuant(s, t, mid, terms);
}

void ContinuantCache::clear() {
	table.clear();
	table_hit.clear();

	available_threads = std::thread::hardware_concurrency() - 1;
}
