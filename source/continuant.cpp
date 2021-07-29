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

mpz_class ContinuantCache::cached_continuant(size_t s, size_t t, size_t mid, const CFTerms& terms) {
	if (table.find({ s,t }) != table.end()) {
		return table[{s, t}];
	}

	auto& ret = table[{s, t}];

	if (t - s <= Params::ThresholdUseBasicContProc)
		ret = basic_continuant(s, t, terms);
	else {

		size_t mid1 = (s + mid) / 2, mid2 = (mid + t) / 2;

		ret = cached_continuant(s, mid, (s + mid) / 2, terms)
			* cached_continuant(mid + 1, t, (mid + t) / 2, terms);
		mpz_addmul(ret.get_mpz_t(),
			cached_continuant(s, mid - 1, (s + mid) / 2, terms).get_mpz_t(),
			cached_continuant(mid + 2, t, (mid + t) / 2, terms).get_mpz_t());
	}

	return ret;
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

mpz_class ContinuantCache::continuant_new(size_t s, size_t t, size_t mid, const CFTerms& terms) {
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
		ret = continuant_new(s, mid, (s + mid) / 2, terms)
			* continuant_new(mid + 1, t, (mid + t) / 2, terms);
		mpz_addmul(ret.get_mpz_t(),
			continuant_new(s, mid - 1, (s + mid) / 2, terms).get_mpz_t(),
			continuant_new(mid + 2, t, (mid + t) / 2, terms).get_mpz_t());
	}

	return ret;
}

mpz_class ContinuantCache::parallel_cached_continuant(size_t s, size_t t, size_t mid, const CFTerms& terms) {
	mu.lock();
	if (table.find({ s,t }) != table.end()) {
		mpz_class ret(table[{s, t}]);
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
			parallel_cached_continuant(s, mid, (s + mid) / 2, terms).get_mpz_t(),
			parallel_cached_continuant(mid + 1, t, (mid + t) / 2, terms).get_mpz_t());
		mpz_addmul(ret.get_mpz_t(),
			parallel_cached_continuant(s, mid - 1, (s + mid) / 2, terms).get_mpz_t(),
			parallel_cached_continuant(mid + 2, t, (mid + t) / 2, terms).get_mpz_t());
	}
	else {

		std::thread th[3];
		bool wait_for_th[3] = {false, false, false};

		if (available_threads > 0) {
			available_threads--;
			wait_for_th[0] = true;

			th[0] = std::thread(parallel_cached_continuant, s, mid, (s + mid) / 2, terms);
		}
		else parallel_cached_continuant(s, mid, (s + mid) / 2, terms);

		if (available_threads > 0) {
			available_threads--;
			wait_for_th[1] = true;

			th[1] = std::thread(parallel_cached_continuant, mid + 1, t, (mid + t) / 2, terms);
		}
		else parallel_cached_continuant(mid + 1, t, (mid + t) / 2, terms);
		
		if (available_threads > 0) {
			available_threads--;
			wait_for_th[2] = true;

			th[2] = std::thread(parallel_cached_continuant, s, mid - 1, (s + mid) / 2, terms);
		}
		else parallel_cached_continuant(s, mid - 1, (s + mid) / 2, terms);

		parallel_cached_continuant(mid + 2, t, (mid + t) / 2, terms);

		for (int i = 0; i < 3; ++i)
			if (wait_for_th[i]) {
				th[i].join();
				available_threads++;
			}

    mu.lock();
    mpz_class u1(table[{s, mid    }]), u2(table[{mid + 1, t}]),
      b1(table[{s, mid - 1}]), b2(table[{mid + 2, t}]);
    mu.unlock();

		mpz_addmul(ret.get_mpz_t(), u1.get_mpz_t(), u2.get_mpz_t());
		mpz_addmul(ret.get_mpz_t(), b1.get_mpz_t(), b2.get_mpz_t());
	}

	mu.lock();
	table[{s, t}] = ret;
	mu.unlock();

	return ret;
}

mpz_class ContinuantCache::parallel_continuant_new(size_t s, size_t t, size_t mid, const CFTerms& terms) {
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
			parallel_continuant_new(s, mid, (s + mid) / 2, terms).get_mpz_t(),
			parallel_continuant_new(mid + 1, t, (mid + t) / 2, terms).get_mpz_t());
		mpz_addmul(ret.get_mpz_t(),
			parallel_continuant_new(s, mid - 1, (s + mid) / 2, terms).get_mpz_t(),
			parallel_continuant_new(mid + 2, t, (mid + t) / 2, terms).get_mpz_t());
	}
	else {
		mpz_class m[4];
		std::thread th;
		auto new_thread_func = [&]() {
			m[1] = parallel_continuant_new(mid + 1, t, (mid + t) / 2, terms);
			m[3] = parallel_continuant_new(mid + 2, t, (mid + t) / 2, terms);
		};

		if (available_threads > 0) {
			available_threads--;
			th = std::thread(new_thread_func);
		}
		else new_thread_func();

		m[0] = parallel_continuant_new(s, mid, (s + mid) / 2, terms);
		m[2] = parallel_continuant_new(s, mid - 1, (s + mid) / 2, terms);

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

	//dummy_run(s, t, mid);
	return parallel_continuant_new(s, t, mid, terms);
	/*if (t - s > Params::ContinuantUseParallel && available_threads >= 3)
		return parallel_continuant(s, t, terms, split_point);

	return cached_continuant(s, mid, mid / 2, terms)
		* cached_continuant(mid + 1, t, (mid + 1 + t) / 2, terms)
		+ cached_continuant(s, mid - 1, mid / 2, terms)
		* cached_continuant(mid + 2, t, (mid + 1 + t) / 2, terms);*/
}

mpz_class parallel_continuant(size_t s, size_t t, const CFTerms& terms, size_t split_point) {
	auto mid = split_point == 0 ? (s + t) / 2 : split_point;

	using namespace ContinuantCache;

	available_threads -= 3;
	std::thread t1(parallel_cached_continuant, s, mid, mid / 2, terms);
	std::thread t2(parallel_cached_continuant, mid + 1, t, (mid + 1 + t) / 2, terms);
	std::thread t3(parallel_cached_continuant, s, mid - 1, mid / 2, terms);
	mpz_class btm(parallel_cached_continuant(mid + 2, t, (mid + 1 + t) / 2, terms));

	t1.join();
	available_threads++;
	t2.join();
	available_threads++;
	t3.join();
	available_threads++;
	
	mpz_class up(table[{s, mid}] * table[{mid + 1, t}]);
	return up + table[{s, mid - 1}] * btm;
}

void ContinuantCache::clear() {
	table.clear();
	table_hit.clear();

	available_threads = std::thread::hardware_concurrency() - 1;
}
