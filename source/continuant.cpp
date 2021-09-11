#include "continuant.h"
#include "mul.h"
#include <cassert>
#include <thread>
#include <mutex>
#include <atomic>
#include <future>

#define USE_YMP_IN_CONTINUANT

namespace ContinuantCache {
	std::map<std::pair<size_t, size_t>, mpz_class> table;
	std::map<std::pair<size_t, size_t>, int> index_hits;

	std::mutex mu;
	std::atomic<int> available_threads;
}

void ContinuantCache::build_cache_indices(size_t s, size_t t, size_t mid) {
	if (index_hits.find({ s,t }) != index_hits.end()) {
		index_hits[{s, t}]++;
		return;
	}

	index_hits[{s, t}] = 0;

	if (t - s <= Params::ThresholdUseBasicContProc)
		return;
	else
	{
		build_cache_indices(s, mid, (s + mid) / 2);
		build_cache_indices(mid + 1, t, (mid + t) / 2);
		build_cache_indices(s, mid - 1, (s + mid) / 2),
		build_cache_indices(mid + 2, t, (mid + t) / 2);
	}
}

mpz_class ContinuantCache::single_threaded_continuant(size_t s, size_t t, size_t mid, const CFTerms& terms)
{
	if (table.find({ s,t }) != table.end())
	{
		index_hits[{s, t}]--;

		// Delete from cache table if no longer usable
		if (index_hits[{s, t}] == 0)
		{
			mpz_class ret(table[{s, t}]);
			table.erase({ s, t });
			return ret;
		}

		return table[{s, t}];
	}

	mpz_class r;
	auto& ret = index_hits[{s, t}] > 0 ? table[{s, t}] : r;

	if (t - s <= Params::ThresholdUseBasicContProc)
		ret = basic_continuant(s, t, terms);
	else {
#ifndef USE_YMP_IN_CONTINUANT
		ret = single_threaded_continuant(s, mid, (s + mid) / 2, terms)
			* single_threaded_continuant(mid + 1, t, (mid + t) / 2, terms);

		mpz_addmul(ret.get_mpz_t(),
			single_threaded_continuant(s, mid - 1, (s + mid) / 2, terms).get_mpz_t(),
			single_threaded_continuant(mid + 2, t, (mid + t) / 2, terms).get_mpz_t());
#else
		multiplication::mul(ret.get_mpz_t(),
			single_threaded_continuant(s, mid, (s + mid) / 2, terms).get_mpz_t(),
			single_threaded_continuant(mid + 1, t, (mid + t) / 2, terms).get_mpz_t());
		multiplication::addmul(ret.get_mpz_t(),
			single_threaded_continuant(s, mid - 1, (s + mid) / 2, terms).get_mpz_t(),
			single_threaded_continuant(mid + 2, t, (mid + t) / 2, terms).get_mpz_t());
#endif
	}

	return ret;
}

mpz_class ContinuantCache::multi_threaded_continuant(size_t s, size_t t, size_t mid, const CFTerms& terms)
{
	mu.lock();
	if (table.find({ s,t }) != table.end())
	{
		mpz_class ret(table[{s, t}]);

		index_hits[{s, t}]--;

		if (index_hits[{s, t}] == 0) {
			table.erase({ s,t });
		}

		mu.unlock();
		return ret;
	}
	mu.unlock();

	mpz_class ret(0);

	if (t - s <= Params::ThresholdUseBasicContProc) ret = basic_continuant(s, t, terms);
	else if (t - s < Params::ContinuantUseParallel || available_threads <= 0)
	{
		mpz_mul(ret.get_mpz_t(),
			multi_threaded_continuant(      s, mid, (s + mid) / 2, terms).get_mpz_t(),
      multi_threaded_continuant(mid + 1,   t, (mid + t) / 2, terms).get_mpz_t()
		);
		mpz_addmul(ret.get_mpz_t(),
			multi_threaded_continuant(      s, mid - 1, (s + mid) / 2, terms).get_mpz_t(),
			multi_threaded_continuant(mid + 2,       t, (mid + t) / 2, terms).get_mpz_t()
		);
	}
	else
	{
		mpz_class m[4];
		std::thread th;
		auto new_thread_func = [&]()
		{
			m[1] = multi_threaded_continuant(mid + 1, t, (mid + t) / 2, terms);
			m[3] = multi_threaded_continuant(mid + 2, t, (mid + t) / 2, terms);
		};

		if (available_threads)
		{
			available_threads--;
			th = std::thread(new_thread_func);
		}
		else new_thread_func();

		m[0] = multi_threaded_continuant(s, mid, (s + mid) / 2, terms);
		m[2] = multi_threaded_continuant(s, mid - 1, (s + mid) / 2, terms);

		if (th.joinable())
		{
			th.join();
			available_threads++;
		}

		mpz_class ret1, ret2;
		if (available_threads)
		{
			available_threads--;
			th = std::thread(mpz_addmul, ret1.get_mpz_t(), m[0].get_mpz_t(), m[1].get_mpz_t());
		}
		else mpz_addmul(ret1.get_mpz_t(), m[0].get_mpz_t(), m[1].get_mpz_t());
		
		mpz_addmul(ret2.get_mpz_t(), m[2].get_mpz_t(), m[3].get_mpz_t());

		if (th.joinable())
		{
			th.join();
			available_threads++;
		}

		mpz_add(ret.get_mpz_t(), ret1.get_mpz_t(), ret2.get_mpz_t());
	}

	mu.lock();
	if (index_hits[{s, t}] > 0) 
		table[{s, t}] = ret;
	mu.unlock();

	return ret;
}

mpz_class continuant(size_t s, size_t t, const CFTerms& terms, size_t split_point)
{
	auto mid = split_point == 0 ? (s+t)/2 : split_point;

	using namespace ContinuantCache;

	if (t - s > Params::ContinuantUseParallel)
		return multi_threaded_continuant (s, t, mid, terms);
	else
		return single_threaded_continuant(s, t, mid, terms);
}

void combine_continuants(CFTermList::Continuants& upperhalf, const CFTermList::Continuants& lowerhalf, bool only_two)
{
	mpz_class p_k, q_k, p_k1, q_k1;

	// TODO Move to Continuant namespace via a combine_continuants() func
	multiplication::mul(p_k.get_mpz_t(), upperhalf.p_k.get_mpz_t(), lowerhalf.p_k.get_mpz_t());
	multiplication::mul(q_k.get_mpz_t(), upperhalf.q_k.get_mpz_t(), lowerhalf.p_k.get_mpz_t());

	if (!only_two) {
		multiplication::mul(p_k1.get_mpz_t(), upperhalf.p_k.get_mpz_t(), lowerhalf.p_k1.get_mpz_t());
		multiplication::mul(q_k1.get_mpz_t(), upperhalf.q_k.get_mpz_t(), lowerhalf.p_k1.get_mpz_t());
	}

	multiplication::addmul(p_k.get_mpz_t(), upperhalf.p_k1.get_mpz_t(), lowerhalf.q_k.get_mpz_t());
	multiplication::addmul(q_k.get_mpz_t(), upperhalf.q_k1.get_mpz_t(), lowerhalf.q_k.get_mpz_t());

	if (!only_two) {
		multiplication::addmul(p_k1.get_mpz_t(), upperhalf.p_k1.get_mpz_t(), lowerhalf.q_k1.get_mpz_t());
		multiplication::addmul(q_k1.get_mpz_t(), upperhalf.q_k1.get_mpz_t(), lowerhalf.q_k1.get_mpz_t());
	}

	upperhalf.p_k = p_k;
	upperhalf.q_k = q_k;
	upperhalf.p_k1 = p_k1;
	upperhalf.q_k1 = q_k1;
}

void ContinuantCache::clear()
{
	table.clear();
	index_hits.clear();

	available_threads = std::thread::hardware_concurrency() - 1;
}
