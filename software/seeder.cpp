#include "graph.h"

#include "tbb/mutex.h"
#include "tbb/parallel_for_each.h"

filter_input seeder_body::operator()(seeder_input input)
{
	reader_output &reads = get<0>(input);

	size_t token = get<1>(input);

	seeder_data output;

	output.fwAnchorBuckets.push_back(0ull);
	output.rcAnchorBuckets.push_back(0ull);

	tbb::mutex output_mutex;

	tbb::parallel_for_each( reads.cbegin(), reads.cend(),
		[&](const Read &read) 
	{
		const size_t read_len = read.seq.size();
		char *read_char = (char *)read.seq.data();
		char *rev_read_char = (char *)read.rc_seq.data();

		uint32_t kmer_max_occurence = cfg.seed_occurence_multiple * (1 + ((g_DRAM->referenceSize) >> (2 * cfg.seed_size)));
		size_t max_hits = kmer_max_occurence * cfg.num_seeds;

		// Forward reads
		auto fwAnchors = sa->DSOFT(read_char, read_len, cfg.num_seeds, cfg.dsoft_threshold, max_hits, cfg.max_candidates, cfg.do_overlap);

		output_mutex.lock();
		output.fwAnchors.insert(output.fwAnchors.end(), fwAnchors.begin(), fwAnchors.end());
		output.fwAnchorBuckets.push_back(output.fwAnchorBuckets.back() + fwAnchors.size());
		output_mutex.unlock();

		// Reverse-complement reads
		auto rcAnchors = sa->DSOFT(rev_read_char, read_len, cfg.num_seeds, cfg.dsoft_threshold, max_hits, cfg.max_candidates, cfg.do_overlap);

		output_mutex.lock();
		output.rcAnchors.insert(output.rcAnchors.end(), rcAnchors.begin(), rcAnchors.end());
		output.rcAnchorBuckets.push_back(output.rcAnchorBuckets.back() + rcAnchors.size());
		output_mutex.unlock();
	});

	return filter_input(filter_payload(reads, output), token);
}
