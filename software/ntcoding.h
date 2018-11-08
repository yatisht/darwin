#pragma once

#define A_NT 0
#define C_NT 1
#define G_NT 2
#define T_NT 3
#define N_NT 4

#include <stdint.h>
#include <iostream>
#include <assert.h>
#include <vector>

#define FIND_MIN(x, y) (x < y) ? x : y
#define FIND_MAX(x, y) (x > y) ? x : y


using namespace std;

void SetIgnoreLower();
uint32_t NtChar2Int (char nt);
void GenerateShapePos(std::string shape);
uint32_t KmerToIndex(std::string kmer);
uint32_t GetKmerIndexAtPos(char* sequence, uint32_t pos);
uint32_t GetKmerIndexAtPos(std::string sequence, uint32_t pos);
static inline int NtToTwoBit (char nt);
uint32_t* SeqToTwoBit (char* seq, uint32_t seq_len);
std::pair<uint64_t*, uint32_t> TwoBitToMinimizers (uint32_t* s_2bit, uint32_t s_len, int k, int w);
std::pair<uint64_t*, uint32_t> QTwoBitToMinimizers (uint32_t* s_2bit, uint32_t s_len, int k, int w);
void SeqToMinimizers(
	std::vector<uint64_t>& minimizers, int64_t &miniCount,
	int kmer_size, int minimizer_window,
	char* seq_addr, std::size_t seq_len, std::size_t seq_offset);

inline uint32_t Min_Window(uint32_t* window, int size) {
	uint32_t min = (1ull << 31ull) - 1ull;
	for (int i = 0; i < size; i++) {
		if (min > window[i]) {
			min = window[i];
		}
	}
	return min;
}

inline uint32_t GetSeedAtPos(uint32_t* s_2bit, uint32_t pos, int k) {
	int idx, shift;
	uint32_t m = (1 << 2 * k) - 1;
	uint32_t seed;
	idx = pos / 16;
	shift = pos % 16;
	uint64_t concat = (((uint64_t)s_2bit[idx + 1]) << 32) + s_2bit[idx];
	seed = (concat >> 2 * shift) & m;
	return seed;
}

inline uint32_t hash32(uint32_t key, int k)
{
	uint32_t m = (1 << 2 * k) - 1;
	key = (~key + (key << 21)) & m;
	key = key ^ (key >> 24);
	key = ((key + (key << 3)) + (key << 8)) & m;
	key = key ^ (key >> 14);
	key = ((key + (key << 2)) + (key << 4)) & m;
	key = key ^ (key >> 28);
	key = (key + (key << 31)) & m;
	return key;
}
