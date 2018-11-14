#define NOMINMAX
#include <algorithm>

#include "tbb/parallel_scan.h"
#include "tbb/parallel_sort.h"
#include "tbb/blocked_range.h"
#include "tbb/scalable_allocator.h"

#include "seed_pos_table.h"
#include "Processor.h"

SeedPosTable::SeedPosTable() {
	ref_size_ = 0;
	kmer_size_ = 0;
	minimizer_window_ = 0;
	max_stride_ = 1;
	shape_size_ = 0;
	num_bins_ = 0;
	avg_hits_ = 0;
}

int SeedPosTable::GetKmerSize() {
	return kmer_size_;
}

bool SeedPosTable::IsPresent(uint32_t index) {
	uint32_t start_index = seedBuckets[index];
	uint32_t end_index = seedBuckets[index + 1];
	return (end_index - start_index <= kmer_max_occurence_);
}


template <size_t __N>
inline __m256i _mm256_shift_left_si256(__m256i a) {
	__m256i c = _mm256_permute2x128_si256(a, a, 0x08);
	return _mm256_alignr_epi8(a, c, 16 - __N);
}

SeedPosTable::SeedPosTable(uint32_t ref_length, const int seed_size, const int minimizer_window, int max_stride, const uint32_t seed_occurence_multiple, const uint32_t bin_size,
	const tbb::concurrent_vector<mini_list> &minimizers,
	uint32_t* seedHistogram, const std::size_t histogramSize) {
	int kmer_size = seed_size;
	shape_size_ = seed_size;
	minimizer_window_ = minimizer_window;
	max_stride_ = max_stride;

	assert(kmer_size <= 15);
	assert(kmer_size > 3);
	kmer_size_ = kmer_size;
	bin_size_ = bin_size;
	log_bin_size_ = (uint32_t)(log2(bin_size_));
	ref_size_ = ref_length;

	kmer_max_occurence_ = seed_occurence_multiple * (1 + ((ref_length) >> (2 * kmer_size)));

	uint32_t seedCount = 0;
	for (auto &miniList : minimizers)
	{
		seedCount += miniList.size();
	}

	// Prefix Sum
	seedBuckets = (uint32_t*)scalable_aligned_malloc((histogramSize + 1) * sizeof(uint32_t), 32);

	__m256i _sum = _mm256_setzero_si256();
	_mm256_storeu_si256((__m256i*)(seedBuckets), _sum);

	for (size_t bucketIndex = 0; bucketIndex < histogramSize; bucketIndex += 8)
	{
		__m256i _buckets = _mm256_load_si256((__m256i*)(seedHistogram + bucketIndex));

		// Vector prefix sum
		_buckets = _mm256_add_epi32(
			_buckets,
			_mm256_shift_left_si256<4>(
				_buckets));

		_buckets = _mm256_add_epi32(
			_buckets,
			_mm256_shift_left_si256<8>(
				_buckets));

		_buckets = _mm256_add_epi32(
			_buckets,
			_mm256_permute2x128_si256(
				_buckets,
				_buckets,
				0x08));

		_buckets = _mm256_add_epi32(
			_buckets,
			_sum);

		_mm256_storeu_si256((__m256i*)(seedBuckets + bucketIndex + 1), _buckets);

		_sum = _mm256_permutevar8x32_epi32(
			_buckets,
			_mm256_setr_epi32(7, 7, 7, 7, 7, 7, 7, 7));
	}

	/*
	uint32_t *seedBuckets2 = (uint32_t*)scalable_malloc((histogramSize + 1) * sizeof(uint32_t));
	uint32_t sum = 0;
	for (size_t i = 0; i < histogramSize; i++)
	{
	seedBuckets2[i] = sum;
	sum += seedHistogram[i];
	}
	seedBuckets2[histogramSize] = sum;

	for (size_t i = 0; i <= histogramSize; i++)
	{
	assert(seedBuckets[i] == seedBuckets2[i]);
	}
	*/

	avg_hits_ = 1 + (seedCount / (1 << 2 * kmer_size_));

	assert(seedBuckets[histogramSize] == seedCount);

	// Bucket sort
	seedPositions = (uint32_t*)scalable_malloc(seedCount * sizeof(uint32_t));

	for (std::size_t chr = 0; chr <  minimizers.size(); chr++)
	{
		auto miniList = minimizers[chr];

		const uint32_t localCount = miniList.size();

		for (uint32_t i = 0; i < localCount; i++)
		{
			uint64_t minimizer = miniList[i];

			uint32_t pos = ((minimizer << 32) >> 32);
			uint32_t seed = (minimizer >> 32);

			uint32_t posIndex = seedBuckets[seed] + (--seedHistogram[seed]);
			seedPositions[posIndex] = pos;
		}
	}

	//tbb::parallel_for(tbb::blocked_range<std::size_t>(0, minimizers.size()),
	//	[&](const tbb::blocked_range<std::size_t> &r) {
	//	for (std::size_t chr = r.begin(); chr < r.end(); chr++)
	//	{
	//		auto miniList = minimizers[chr];

	//		const uint32_t localCount = miniList.size();

	//		for (uint32_t i = 0; i < localCount; i++)
	//		{
	//			uint64_t minimizer = miniList[i];

	//			uint32_t pos = ((minimizer << 32) >> 32);
	//			uint32_t seed = (minimizer >> 32);

	//			uint32_t posIndex = seedBuckets[seed] + InterlockedDecrement(seedHistogram + seed);
	//			seedPositions[posIndex] = pos;
	//		}
	//	}
	//});

	// sort buckets by position
	tbb::parallel_for(tbb::blocked_range<std::size_t>(0, histogramSize),
		[&](const tbb::blocked_range<std::size_t> &r) {
		for (std::size_t i = r.begin(); i < r.end(); i++) {
			if ((seedBuckets[i] < seedBuckets[i + 1]) && (seedBuckets[i + 1] - seedBuckets[i] <= kmer_max_occurence_))
			{
				std::sort(seedPositions + seedBuckets[i], seedPositions + seedBuckets[i + 1]);
			}
		}
	});

	/*
	// Histogram should be empty by now
	for (uint32_t i = 0; i < histogramSize; i++) {
	assert(seedHistogram[i] == 0);
	}

	uint64_t* miniArray = (uint64_t*)scalable_malloc(ref_length * sizeof(uint64_t));

	// Compact the minimizers
	uint32_t curr = 0;

	for (auto &miniList : minimizers)
	{
	const uint32_t localCount = miniList.size();

	for (uint32_t i = 0; i < localCount; i++)
	{
	miniArray[curr++] = miniList[i];
	}
	}

	assert(curr == seedCount);

	// Compare with the slow sort method
	tbb::parallel_sort(miniArray, miniArray + seedCount);

	uint32_t *index_table_;
	uint32_t *pos_table_;

	index_table_size_ = ((uint32_t)1 << 2 * kmer_size) + 1;
	index_table_ = new uint32_t[index_table_size_];

	pos_table_ = new uint32_t[seedCount];

	uint32_t curr_index = 0;
	uint32_t seed, pos;

	for (uint32_t i = 0; i < seedCount; i++) {
	pos = ((miniArray[i] << 32) >> 32);
	seed = (miniArray[i] >> 32);
	pos_table_[i] = pos;
	if (seed > curr_index) {
	for (uint32_t s = curr_index; s < seed; s++) {
	index_table_[s] = i;
	}
	curr_index = seed;
	}
	}
	for (uint32_t i = curr_index; i < index_table_size_; i++) {
	index_table_[i] = seedCount;
	}

	for (uint32_t i = 0; i < histogramSize; i++) {
	assert(index_table_[i] == seedBuckets[i + 1]);
	}

	for (uint32_t i = 0; i < seedCount; i++) {
	assert(pos_table_[i] == seedPositions[i]);
	}

	delete[] index_table_;
	delete[] pos_table_;
	scalable_free(miniArray);
	*/
}

SeedPosTable::~SeedPosTable() {
	scalable_aligned_free(seedBuckets);
	scalable_free(seedPositions);
}

static inline int collinear_score(uint32_t h1, uint32_t o1, uint32_t h2, uint32_t o2) {
	return ((h1 <= h2) && (o1 <= o2)) ? 1 : -(1 << 24);
}


std::vector<Anchors> SeedPosTable::DSOFT(char* query, uint32_t query_length, int N, int threshold, size_t max_hits, size_t max_candidates, int overlap)
{
	uint64_t k = kmer_size_;
	uint64_t w = minimizer_window_;

	uint64_t* minimizers = (uint64_t*)scalable_malloc(query_length * sizeof(uint64_t));

	uint64_t min_count = 0;

	iterate_minimizers(query, query_length, k, w,
		[&](uint64_t p, uint32_t m) {
		minimizers[min_count] = (p << 32) + m;
		min_count++;
	});

	//#ifdef _DEBUG
	//	uint32_t qlen_2bit = (query_length + 15) / 16;
	//	uint32_t* q_2bit = SeqToTwoBit(query, query_length);
	//
	//	uint32_t* window = (uint32_t*)calloc(w, sizeof(uint32_t));
	//	uint64_t last_m = 0;
	//	uint32_t last_p = 0;
	//
	//	uint32_t num_min = 0;
	//
	//	for (int p = 0; p < w - 1; p++) {
	//		window[p] = hash32(GetSeedAtPos(q_2bit, p, k), k);
	//	}
	//
	//	for (uint32_t p = w - 1; p < 16 * qlen_2bit - k; p++) {
	//		window[p%w] = hash32(GetSeedAtPos(q_2bit, p, k), k);
	//		uint64_t m = Min_Window(window, w);
	//		if ((m != last_m) || (p - last_p >= w)) {
	//			assert(minimizers[num_min] == ((uint64_t)p << 32) + m);
	//			num_min++;
	//			last_m = m;
	//			last_p = p;
	//		}
	//	}
	//
	//	assert(num_min == min_count);
	//
	//	free(window);
	//#endif

	uint32_t expected_hits = avg_hits_ * min_count;

	std::vector<Hits> hits_array;
	hits_array.reserve(expected_hits);

	int num_seeds = 0;

	int stride = 1; //FIND_MAX(1, FIND_MIN(min_count/N, max_stride_));

	for (int i = 0; i < min_count; i += stride) {

		uint32_t offset = (minimizers[i] >> 32);
		uint32_t index = ((minimizers[i] << 32) >> 32);

		uint32_t start_index = seedBuckets[index];
		uint32_t end_index = seedBuckets[index + 1];

		if (end_index - start_index <= kmer_max_occurence_) {
			num_seeds++;
			for (uint32_t j = start_index; j < end_index; j++) {
				uint32_t hit = seedPositions[j];
				if (hit >= offset) {
					uint32_t bin = ((hit - offset) / bin_size_);
					uint64_t bin_offset = ((uint64_t)bin << 32) + offset;
					hits_array.push_back(Hits(bin_offset, hit));
				}
			}
		}


		if (i > N) {
            if (overlap == 0) {
                stride = max_stride_;
            }
            else {
                break;
            }
		}
	}

	std::stable_sort(hits_array.begin(), hits_array.end(), CompareHits);


	int num_hits = 0;
	int num_candidates = 0;

	std::vector<Anchors> anchors;
	anchors.reserve(max_candidates);
	std::vector<Anchors> final_anchors;
	final_anchors.reserve(max_candidates);
	std::deque<uint32_t> candidate_bins;
	anchors.reserve(max_candidates);

	num_hits = hits_array.size();

	uint32_t last_bin = (1 << 31);
	uint32_t last_offset = 0;
	uint32_t curr_count = 0;

	for (int i = 0; i < num_hits; i++) {
		Hits next_hit = hits_array[i];
		uint32_t offset = ((next_hit.bin_offset << 32) >> 32);
		uint32_t bin = (next_hit.bin_offset >> 32);
		uint32_t hit = next_hit.hit;
		if (bin == last_bin) {
			if (curr_count < threshold) {
				curr_count = ((offset - last_offset > kmer_size_) || (curr_count == 0)) ? curr_count + kmer_size_ : curr_count + (offset - last_offset);
				//curr_count += 1; 
				if (curr_count >= threshold) {
					uint64_t anchor = ((uint64_t)hit << 32) + offset;
					anchors.push_back(Anchors(anchor));
					candidate_bins.push_back(bin);
					if (num_candidates >= max_candidates) {
						break;
					}
				}
			}
		}
		else {
			last_bin = bin;
			curr_count = kmer_size_;
			//curr_count = 1;
		}
		last_offset = offset;
	}

	int chain_len = 0;
	int max_chain_len = 0;
	int curr_candidate = 0;
	uint32_t sv_num_bins = (overlap == 0) ? (1 << 12) / bin_size_ : 1;

	int start_idx = 0;

	// Each cadidate bin gets a set of hit offsets within a SV window
	for (int k = 0; k < anchors.size(); k++) {
		uint32_t curr_bin = candidate_bins[k];
		bool start_assigned = false;
		anchors[k].num_chained_hits = 0;
		anchors[k].anchor_score = 0;

		for (int i = start_idx; i < num_hits; i++) {
			Hits next_hit = hits_array[i];
			uint32_t bin = (next_hit.bin_offset >> 32);

			if ((bin + sv_num_bins >= curr_bin) && (bin < curr_bin + sv_num_bins)) {
				if (!start_assigned) {
					start_assigned = true;
					start_idx = i;
				}
				uint32_t offset = ((next_hit.bin_offset << 32) >> 32);
				uint32_t hit = next_hit.hit;

				uint64_t hit_offset = ((uint64_t)hit << 32) + offset;
				if (hit_offset <= anchors[k].hit_offset) {
					anchors[k].left_chained_hits.push_back(hit_offset);
				}
				if (hit_offset >= anchors[k].hit_offset) {
					anchors[k].right_chained_hits.push_back(hit_offset);
				}
			}
			else if ((bin >= curr_bin + sv_num_bins)) {
				break;
			}
		}

		assert(anchors[k].left_chained_hits.size() > 0);
		assert(anchors[k].right_chained_hits.size() > 0);
		std::sort(anchors[k].left_chained_hits.begin(), anchors[k].left_chained_hits.end());
		std::sort(anchors[k].right_chained_hits.begin(), anchors[k].right_chained_hits.end());

		// Left chain collinear
		int left_chain_size = anchors[k].left_chained_hits.size();
		vector<uint64_t> collinear_left_chained_hits;
		collinear_left_chained_hits.reserve(left_chain_size);

		uint64_t curr_hit_offset = anchors[k].left_chained_hits.back();
		collinear_left_chained_hits.push_back(curr_hit_offset);

		for (int h = left_chain_size - 2; h >= 0; h--) {
			uint64_t hit_offset = anchors[k].left_chained_hits[h];

			uint32_t h1 = (curr_hit_offset >> 32);
			uint32_t o1 = ((curr_hit_offset << 32) >> 32);
			uint32_t h2 = (hit_offset >> 32);
			uint32_t o2 = ((hit_offset << 32) >> 32);

			if ((h1 >= h2) && (o1 >= o2)) {
				//int match = FIND_MIN((h1 - h2), (o1 - o2));
				int match = ((h1 - h2) > (o1 - o2)) ? (o1 - o2) : (h1 - h2);
				int gap = ((h1 - h2) > (o1 - o2)) ? ((h1 - h2) - (o1 - o2)) : ((o1 - o2) - (h1 - h2));
				anchors[k].anchor_score += (match - gap / 10);
				collinear_left_chained_hits.push_back(hit_offset);
				curr_hit_offset = hit_offset;
			}
		}

		std::sort(collinear_left_chained_hits.begin(), collinear_left_chained_hits.end());
		anchors[k].left_chained_hits.clear();
		anchors[k].left_chained_hits.assign(collinear_left_chained_hits.begin(), collinear_left_chained_hits.end());

		// Right chain collinear
		int right_chain_size = anchors[k].right_chained_hits.size();
		vector<uint64_t> collinear_right_chained_hits;
		collinear_right_chained_hits.reserve(right_chain_size);

		curr_hit_offset = anchors[k].right_chained_hits.front();
		collinear_right_chained_hits.push_back(curr_hit_offset);
		for (int h = 1; h < right_chain_size; h++) {
			uint64_t hit_offset = anchors[k].right_chained_hits[h];

			uint32_t h1 = (curr_hit_offset >> 32);
			uint32_t o1 = ((curr_hit_offset << 32) >> 32);
			uint32_t h2 = (hit_offset >> 32);
			uint32_t o2 = ((hit_offset << 32) >> 32);

			if ((h1 <= h2) && (o1 <= o2)) {
				int match = ((h2 - h1) > (o2 - o1)) ? (o2 - o1) : (h2 - h1);
				int gap = ((h2 - h1) > (o2 - o1)) ? ((h2 - h1) - (o2 - o1)) : ((o2 - o1) - (h2 - h1));
				anchors[k].anchor_score += (match - gap / 10);
				collinear_right_chained_hits.push_back(hit_offset);
				curr_hit_offset = hit_offset;
			}
		}
		std::reverse(collinear_right_chained_hits.begin(), collinear_right_chained_hits.end());
		anchors[k].right_chained_hits.clear();
		anchors[k].right_chained_hits.assign(collinear_right_chained_hits.begin(), collinear_right_chained_hits.end());


		chain_len = anchors[k].left_chained_hits.size() + anchors[k].right_chained_hits.size();
		anchors[k].num_chained_hits = chain_len;
		if (chain_len > max_chain_len) {
			max_chain_len = chain_len;
		}
	}


	int num_secondary_chains = 0;

	// Discard anchors with overlapping range early
	if (!anchors.empty())
	{
		std::sort(anchors.begin(), anchors.end(),
			[](const Anchors &a, const Anchors &b) {
			//return ((a.anchor_score > b.anchor_score) || ((a.anchor_score == b.anchor_score) && (a.hit_offset < b.hit_offset)));
			return ((a.num_chained_hits > b.num_chained_hits) || ((a.num_chained_hits == b.num_chained_hits) && (a.hit_offset < b.hit_offset)));
		});

		for (auto a1 = anchors.begin(); a1 != anchors.end(); a1++)
		{
			uint32_t hs1 = (a1->left_chained_hits.front() >> 32);
			uint32_t he1 = (a1->right_chained_hits.front() >> 32);
			uint32_t s1 = ((a1->left_chained_hits.front() << 32) >> 32);
			uint32_t e1 = ((a1->right_chained_hits.front() << 32) >> 32);

			if (a1->num_chained_hits > 0) {
				num_secondary_chains++;
			}

			a1->num_chained_hits = a1->left_chained_hits.size() + a1->right_chained_hits.size();

			final_anchors.push_back(*a1);

			if ((num_secondary_chains > 5) && (overlap == 1)) {
				break;
			}

            if (overlap == 1) {
                for (auto a2 = a1 + 1; a2 != anchors.end(); a2++)
                {
                    float r1 = (a1->hit_offset >> 32);
                    float q1 = ((a1->hit_offset << 32) >> 32);
                    float r2 = (a2->hit_offset >> 32);
                    float q2 = ((a2->hit_offset << 32) >> 32);

                    float slope = std::abs((r1 - r2) / (q1 - q2) - 1);
                    if (std::abs((r1 - r2) / (q1 - q2) - 1) <= 0.15)
                    {
                        a2->num_chained_hits = -1;
                    }
                }
            }
		}

	}

	scalable_free(minimizers);

	return final_anchors;
}

