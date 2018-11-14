#pragma once

#include <string>
#include <vector>
#include <map>
#include <assert.h>
#include <stdint.h>
#include <math.h>

#include "immintrin.h"

#include "tbb/concurrent_vector.h"

#include "ntcoding.h"

using namespace std;

#define nz_bins 25000000

struct Hits {
    Hits(uint64_t a, uint32_t b)
        : bin_offset(a),
        hit(b)
    {};

    uint64_t bin_offset;
    uint32_t hit;
};

struct Anchors {
    Anchors(uint64_t a)
        : hit_offset(a)
    {};

    uint64_t hit_offset;
    int num_chained_hits;
    int anchor_score;
    vector<uint64_t> left_chained_hits;
	vector<uint64_t> right_chained_hits;
};

static inline bool CompareHits (Hits h1, Hits h2) {
    //return ((h1.bin_offset < h2.bin_offset) || ((h1.bin_offset == h2.bin_offset) && (h1.hit < h2.hit)));
	return ((h1.bin_offset < h2.bin_offset));
}

typedef std::vector<uint64_t> mini_list;

#define _mm256_shuffle_epx32(a, b, imm) _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(a), _mm256_castsi256_ps(b), imm))

template <size_t __N>
inline __m256i _mm256_shift_left_si256(__m256i a, __m256i b) {
	__m256i c = _mm256_permute2x128_si256(a, b, 0x03);
	return _mm256_alignr_epi8(a, c, 16 - __N);
}

template <size_t __N>
inline __m256i _mm256_shift_right_si256(__m256i a, __m256i b) {
	__m256i c = _mm256_permute2x128_si256(a, b, 0x21);
	return _mm256_alignr_epi8(c, a, __N);
}

inline uint64_t _q_to_2_bit(const char* pQuery)
{
	// Load 16 characters 
	__m128i _q = _mm_loadu_si128((const __m128i*) pQuery);

	__m128i _q_2bit = _mm_shuffle_epi8(
		// 8 -> 2 bit conversion table, N character is mapped to zero
		_mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 3, 1, 0, 0, 0),
		_mm_and_si128(
			_q,
			// 8 -> 2 bit conversion mask
            _mm_set1_epi8(0x0f)));

	// PEXT mask for compacting the sequence
	const uint64_t pext_mask_8_2 = 0x0303030303030303ull;

	uint64_t q_2bit =
		_pext_u64(_mm_extract_epi64(_q_2bit, 1), pext_mask_8_2) << 16 |
		_pext_u64(_mm_extract_epi64(_q_2bit, 0), pext_mask_8_2);

	return q_2bit;
}

inline __m256i _get_seeds(const __m256i _min_qword, const __m256i _kmer_mask)
{
	// Vector of 8 minimizers
	__m256i _min_x_8 = _mm256_and_si256(
        _mm256_shuffle_epx32(
			_mm256_srlv_epi64(
				_min_qword,
				_mm256_set_epi64x(
					10ull, 8ull, 2ull, 0ull)),
			_mm256_srlv_epi64(
				_min_qword,
				_mm256_set_epi64x(
					14ull, 12ull, 6ull, 4ull)),
			0x88),
		_kmer_mask);

	// hash32
	//key = (~key + (key << 21)) & m;
	_min_x_8 = _mm256_and_si256(
		_mm256_add_epi32(
			_mm256_xor_si256(
				_min_x_8,
				_mm256_set1_epi32(-1)),
			_mm256_slli_epi32(
				_min_x_8,
				21)),
		_kmer_mask);

	//key = key ^ (key >> 24);
	_min_x_8 = _mm256_xor_si256(
		_min_x_8,
		_mm256_srli_epi32(
			_min_x_8,
			24));

	//key = ((key + (key << 3)) + (key << 8)) & m;
	_min_x_8 = _mm256_and_si256(
		_mm256_add_epi32(
			_min_x_8,
			_mm256_add_epi32(
				_mm256_slli_epi32(
					_min_x_8,
					3),
				_mm256_slli_epi32(
					_min_x_8,
					8))),
		_kmer_mask);

	//key = key ^ (key >> 14);
	_min_x_8 = _mm256_xor_si256(
		_min_x_8,
		_mm256_srli_epi32(
			_min_x_8,
			14));

	//key = ((key + (key << 2)) + (key << 4)) & m;
	_min_x_8 = _mm256_and_si256(
		_mm256_add_epi32(
			_min_x_8,
			_mm256_add_epi32(
				_mm256_slli_epi32(
					_min_x_8,
					2),
				_mm256_slli_epi32(
					_min_x_8,
					4))),
		_kmer_mask);

	//key = key ^ (key >> 28);
	_min_x_8 = _mm256_xor_si256(
		_min_x_8,
		_mm256_srli_epi32(
			_min_x_8,
			28));

	//key = (key + (key << 31)) & m;
	_min_x_8 = _mm256_and_si256(
		_mm256_add_epi32(
			_min_x_8,
			_mm256_slli_epi32(
				_min_x_8,
				31)),
		_kmer_mask);

	return _min_x_8;
}

inline __m256i _sliding_minimum_5w(__m256i &_min_window, const __m256i _min_seeds)
{
	static const __m256i _mm256_000f_epu32 = _mm256_set_epi32(
		0, 0, 0, 0xffffffff,
		0, 0, 0, 0xffffffff);

	// Suffix scan
	__m256i _min_suffix = _mm256_min_epu32(
		_min_seeds,
		_mm256_or_si256(
			_mm256_bsrli_epi128(
				_min_seeds,
				4),
			_mm256_bslli_epi128(
				_mm256_000f_epu32,
				12)));

	_min_suffix = _mm256_min_epu32(
		_min_suffix,
        _mm256_shuffle_epx32(
			_min_suffix,
			_mm256_000f_epu32,
			0x0e));

	// Prefix scan
	__m256i _min_prefix = _mm256_min_epu32(
		_mm256_or_si256(
			_mm256_bslli_epi128(
				_min_seeds,
				4),
			_mm256_000f_epu32),
		_min_seeds);

	_min_prefix = _mm256_min_epu32(
        _mm256_shuffle_epx32(
			_mm256_000f_epu32,
			_min_prefix,
			0x40),
		_min_prefix);

	// Sliding window minimum
	__m256i _min_result = _mm256_min_epu32(
		_mm256_permute2x128_si256(
			_min_window,
			_min_suffix,
			0x21),
		_min_prefix);

	_min_window = _min_suffix;

	return _min_result;
}

inline __m256i _sliding_minimum_9w(__m256i &_min_window, const __m256i _min_seeds)
{
	static const __m256i _mm256_ff_epu32 = _mm256_set1_epi32(0xffffffff);

	// Suffix scan
	__m256i _min_suffix = _mm256_min_epu32(
		_min_seeds,
		_mm256_shift_right_si256<4>(
			_min_seeds,
			_mm256_ff_epu32));

	_min_suffix = _mm256_min_epu32(
		_min_suffix,
		_mm256_shift_right_si256<8>(
			_min_suffix,
			_mm256_ff_epu32));

	_min_suffix = _mm256_min_epu32(
		_min_suffix,
		_mm256_permute2x128_si256(
			_min_suffix,
			_mm256_ff_epu32,
			0x21));

	// Prefix scan
	__m256i _min_prefix = _mm256_min_epu32(
		_min_seeds,
		_mm256_shift_left_si256<4>(
			_min_seeds,
			_mm256_ff_epu32));

	_min_prefix = _mm256_min_epu32(
		_min_prefix,
		_mm256_shift_left_si256<8>(
			_min_prefix,
			_mm256_ff_epu32));

	_min_prefix = _mm256_min_epu32(
		_min_prefix,
		_mm256_permute2x128_si256(
			_min_prefix,
			_mm256_ff_epu32,
			0x03));

	// Sliding window minimum
	__m256i _min_result = _mm256_min_epu32(
		_min_window,
		_min_prefix);

	_min_window = _min_suffix;

	return _min_result;
}

template<typename Body>
void iterate_minimizers_qw(const char* query, const uint32_t query_length, uint64_t k, uint64_t w, Body body)
{
	uint64_t kmer_mask = (1ull << (k << 1)) - 1ull;

	uint64_t last_m = 0, last_p = 0;

	uint32_t* pWindow = (uint32_t*)calloc(w, sizeof(uint32_t));

	// Priming the vector loop
	uint64_t min_current = _q_to_2_bit(query);

	uint64_t min_next = _q_to_2_bit(query + 16);

	uint64_t min_qword = min_current | (min_next << 32);

	uint64_t p = 0;

	uint32_t qlen_centinel = (~0x0f & (query_length + 15)) - k;

	// Fill the window
	for (; p < w - 1; p++) {
		pWindow[p] = hash32(min_qword & kmer_mask, k);
		min_qword >>= 2;
	}

	// Finish the first batch
	for (/*p = w - 1*/; p < 16; p++)
	{
		pWindow[p%w] = hash32(min_qword & kmer_mask, k);
		min_qword >>= 2;
		uint32_t m = Min_Window(pWindow, w);
		if ((m != last_m) || (p - last_p >= w))
		{
			body(p, m);
			last_m = m;
			last_p = p;
		}
	}

	min_current = min_next;

	// Main loop
	while (p < qlen_centinel)
	{
		min_next = _q_to_2_bit(query + p + 16);

		min_qword = min_current | (min_next << 32);

		__m256i _kmer_mask = _mm256_set1_epi32(
			(1 << (k << 1)) - 1);

		__m256i _min_qword = _mm256_set1_epi64x(min_qword);

		__m256i _min_seeds = _get_seeds(_min_qword, _kmer_mask);

		for (uint32_t i = 0; (i < 8) && (p < qlen_centinel); i++, p++)
		{
            __m128i _i = _mm_cvtsi32_si128(i);
            __m256i _val = _mm256_permutevar8x32_epi32(_min_seeds, _mm256_castsi128_si256(_i));
            pWindow[p%w] = _mm_cvtsi128_si32(_mm256_castsi256_si128(_val));

            uint32_t m = Min_Window(pWindow, w);
			if ((m != last_m) || (p - last_p >= w))
			{
				body(p, m);
				last_m = m;
				last_p = p;
			}
		}

		_min_qword = _mm256_srli_epi64(
			_min_qword,
			16);

		_min_seeds = _get_seeds(_min_qword, _kmer_mask);

		for (uint32_t i = 0; (i < 8) && (p < qlen_centinel); i++, p++)
		{
            __m128i _i = _mm_cvtsi32_si128(i);
            __m256i _val = _mm256_permutevar8x32_epi32(_min_seeds, _mm256_castsi128_si256(_i));
            pWindow[p%w] = _mm_cvtsi128_si32(_mm256_castsi256_si128(_val));

			uint32_t m = Min_Window(pWindow, w);
			if ((m != last_m) || (p - last_p >= w))
			{
				body(p, m);
				last_m = m;
				last_p = p;
			}
		}

		min_current = min_next;
	}
}

// Special case with w set to 5
template<typename Body>
void iterate_minimizers_5w(const char* query, const uint32_t query_length, uint64_t k, Body body)
{
	const int64_t w = 5ull;

	uint64_t last_m = 0, last_p = 0;

	// Priming the vector loop
	__m256i _kmer_mask = _mm256_set1_epi32(
		(1 << (k << 1)) - 1);

	uint64_t min_current = _q_to_2_bit(query);

	__m256i _min_window = _mm256_setzero_si256();

	uint64_t qlen_centinel = (~0x0full & (query_length + 15ull)) - k;

	// Main loop
	for (uint64_t p = 0; p < qlen_centinel;)
	{
		uint64_t min_next = _q_to_2_bit(query + p + 16);

		uint64_t min_qword = min_current | (min_next << 32);

		__m256i _min_qword = _mm256_set1_epi64x(min_qword);

		__m256i _min_seeds = _get_seeds(_min_qword, _kmer_mask);

		__m256i _min_result = _sliding_minimum_5w(_min_window, _min_seeds);

		for (uint32_t i = 0; (i < 8) && (p < qlen_centinel); i++, p++)
		{
            __m128i _i = _mm_cvtsi32_si128(i);
            __m256i _val = _mm256_permutevar8x32_epi32(_min_result, _mm256_castsi128_si256(_i));
            uint32_t m = _mm_cvtsi128_si32(_mm256_castsi256_si128(_val));

			if ((m != last_m) || (p - last_p >= w))
			{
				body(p, m);
				last_m = m;
				last_p = p;
			}
		}

		_min_qword = _mm256_srli_epi64(
			_min_qword,
			16);

		_min_seeds = _get_seeds(_min_qword, _kmer_mask);

		_min_result = _sliding_minimum_5w(_min_window, _min_seeds);

		for (uint32_t i = 0; (i < 8) && (p < qlen_centinel); i++, p++)
		{
            __m128i _i = _mm_cvtsi32_si128(i);
            __m256i _val = _mm256_permutevar8x32_epi32(_min_result, _mm256_castsi128_si256(_i));
            uint32_t m = _mm_cvtsi128_si32(_mm256_castsi256_si128(_val));

			if ((m != last_m) || (p - last_p >= w))
			{
				body(p, m);
				last_m = m;
				last_p = p;
			}
		}

		min_current = min_next;
	}
}

// Special case with w set to 9
template<typename Body>
void iterate_minimizers_9w(const char* query, const uint32_t query_length, uint64_t k, Body body)
{
	const int64_t w = 9ull;

	uint64_t last_m = 0, last_p = 0;

	// Priming the vector loop
	__m256i _kmer_mask = _mm256_set1_epi32(
		(1 << (k << 1)) - 1);

	uint64_t min_current = _q_to_2_bit(query);

	__m256i _min_window = _mm256_setzero_si256();

	uint64_t qlen_centinel = (~0x0full & (query_length + 15ull)) - k;

	// Main loop
	for (uint64_t p = 0; p < qlen_centinel;)
	{
		uint64_t min_next = _q_to_2_bit(query + p + 16);

		uint64_t min_qword = min_current | (min_next << 32);

		__m256i _min_qword = _mm256_set1_epi64x(min_qword);

		__m256i _min_seeds = _get_seeds(_min_qword, _kmer_mask);

		__m256i _min_result = _sliding_minimum_9w(_min_window, _min_seeds);

		for (uint32_t i = 0; (i < 8) && (p < qlen_centinel); i++, p++)
		{
            __m128i _i = _mm_cvtsi32_si128(i);
            __m256i _val = _mm256_permutevar8x32_epi32(_min_result, _mm256_castsi128_si256(_i));
            uint32_t m = _mm_cvtsi128_si32(_mm256_castsi256_si128(_val));

			if ((m != last_m) || (p - last_p >= w))
			{
				body(p, m);
				last_m = m;
				last_p = p;
			}
		}

		_min_qword = _mm256_srli_epi64(
			_min_qword,
			16);

		_min_seeds = _get_seeds(_min_qword, _kmer_mask);

		_min_result = _sliding_minimum_9w(_min_window, _min_seeds);

		for (uint32_t i = 0; (i < 8) && (p < qlen_centinel); i++, p++)
		{
            __m128i _i = _mm_cvtsi32_si128(i);
            __m256i _val = _mm256_permutevar8x32_epi32(_min_result, _mm256_castsi128_si256(_i));
            uint32_t m = _mm_cvtsi128_si32(_mm256_castsi256_si128(_val));

			if ((m != last_m) || (p - last_p >= w))
			{
				body(p, m);
				last_m = m;
				last_p = p;
			}
		}

		min_current = min_next;
	}
}

template<typename Body>
void iterate_minimizers(const char* query, const uint32_t query_length, uint64_t k, uint64_t w, Body body)
{
	switch (w)
	{
	case 5:
		iterate_minimizers_5w(query, query_length, k, body);
		break;

	case 9:
		iterate_minimizers_9w(query, query_length, k, body);
		break;

	default:
		iterate_minimizers_qw(query, query_length, k, w, body);
		break;
	}
}

class SeedPosTable {
    private:
        uint32_t index_table_size_;
        uint32_t ref_size_;
        uint32_t bin_size_;
        uint32_t log_bin_size_;
        int kmer_size_;
        int minimizer_window_;
        int max_stride_;
        int shape_size_;
        uint32_t kmer_max_occurence_;

        uint32_t *seedBuckets;
        uint32_t *seedPositions;

        uint32_t num_bins_;

        uint32_t avg_hits_;


    public:
        SeedPosTable();
        SeedPosTable(uint32_t ref_length, const int seed_size, const int minimizer_window, int max_stride, 
                const uint32_t seed_occurence_multiple, const uint32_t bin_size,
                const tbb::concurrent_vector<mini_list> &minimizers,
                uint32_t* seedHistogram, const std::size_t histogramSize);
        //SeedPosTable(char* ref_str, uint32_t ref_length, int seed_size, int minimizer_window, uint32_t seed_occurence_multiple, uint32_t bin_size);
        ~SeedPosTable();

        bool IsPresent(uint32_t index);
        int GetKmerSize();
        std::vector<Anchors> DSOFT(char* query, uint32_t query_length, int N, int threshold, size_t max_hits, size_t max_candidates, int overlap);
};



