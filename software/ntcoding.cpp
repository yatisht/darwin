
#include "ntcoding.h"
#include <cstring>

#include "tbb/scalable_allocator.h"

int shape_pos[32];
int shape_size;
bool is_ignore_lower = false;

uint32_t NtChar2Int(char nt) {
	if (!is_ignore_lower) {
		nt = toupper(nt);
	}
	switch (nt) {
	case 'A': return A_NT;
	case 'C': return C_NT;
	case 'G': return G_NT;
	case 'T': return T_NT;
	case 'N': return N_NT;
	default: return N_NT;
	}
}

void SetIgnoreLower() {
	is_ignore_lower = true;
}

uint32_t KmerToIndex(std::string kmer) {
	uint32_t index = 0;
	char nt;
	for (int i = 0; i < kmer.length(); i++) {
		nt = (kmer.at(i));
		index = (index << 2) + NtChar2Int(nt);
	}
	return index;
}

void GenerateShapePos(std::string shape) {
	shape_size = 0;
	for (int i = 0; i < shape.length(); i++) {
		if (shape[i] == '1') {
			shape_pos[shape_size++] = i;
		}
	}
}

uint32_t GetKmerIndexAtPos(char* sequence, uint32_t pos) {
	uint32_t kmer = 0;
	for (int i = 0; i < shape_size; i++) {
		uint32_t nt = NtChar2Int(sequence[pos + shape_pos[i]]);
		if (nt != N_NT) {
			kmer = (kmer << 2) + nt;
		}
		else {
			kmer = (1 << 31);
			break;
		}
	}
	return kmer;
}

uint32_t GetKmerIndexAtPos(std::string sequence, uint32_t pos) {
	uint32_t kmer = 0;
	for (int i = 0; i < shape_size; i++) {
		uint32_t nt = NtChar2Int(sequence[pos + shape_pos[i]]);
		if (nt != N_NT) {
			kmer = (kmer << 2) + nt;
		}
		else {
			kmer = (1 << 31);
			break;
		}
	}
	return kmer;
}


static inline int NtToTwoBit(char nt) {
	switch (nt) {
	case 'a':
	case 'A': return 0;
	case 'c':
	case 'C': return 1;
	case 'g':
	case 'G': return 2;
	case 't':
	case 'T': return 3;
	default: return 0;
	}
	return 0;
}

uint32_t* SeqToTwoBit(char* seq, uint32_t seq_len) {
	uint32_t len_2bit = (seq_len + 15) / 16;
	uint32_t* s_2bit = (uint32_t*)scalable_calloc(len_2bit, sizeof(uint32_t));

	char c;
	int idx, j_max;
	for (uint32_t i = 0; i < seq_len; i += 16) {
		j_max = FIND_MIN(16, seq_len - i);
		idx = i / 16;
		for (int j = 0; j < j_max; j++) {
			c = seq[i + j];
			s_2bit[idx] += (NtToTwoBit(c) << 2 * j);
		}
	}
	return s_2bit;

}


void SeqToMinimizers(
	std::vector<uint64_t>& minimizers, int64_t &miniCount,
	int kmer_size, int minimizer_window,
	char* seq_addr, std::size_t seq_len, std::size_t seq_offset) {
	uint32_t* window = (uint32_t*)scalable_calloc(minimizer_window, sizeof(uint32_t));
	uint64_t last_m = 0;
	uint32_t last_p = 0;

	std::size_t rlen_2bit = 1 + ((seq_len + 15) / 16);
	uint32_t* r_2bit = SeqToTwoBit(seq_addr, seq_len);

	for (int p = 0; p < minimizer_window - 1; p++) {
		window[p] = hash32(GetSeedAtPos(r_2bit, p, kmer_size), kmer_size);
	}

	for (uint32_t p = minimizer_window - 1; p < seq_len - kmer_size - minimizer_window; p++) {
		window[p%minimizer_window] = hash32(GetSeedAtPos(r_2bit, p, kmer_size), kmer_size);
		uint64_t m = Min_Window(window, minimizer_window);
		if ((m != last_m) || (p - last_p >= minimizer_window))
		{
			int64_t miniPos = ++miniCount;

			minimizers[miniPos] = (m << 32) + p + seq_offset;
			last_m = m;
			last_p = p;
		}
	}

	scalable_free(r_2bit);
	scalable_free(window);
}

std::pair<uint64_t*, uint32_t> TwoBitToMinimizers(uint32_t* s_2bit, uint32_t s_len, int k, int w) {
	uint32_t* window = (uint32_t*)calloc(w, sizeof(uint32_t));
	uint64_t last_m = 0;
	uint32_t last_p = 0;
	uint64_t m;

	uint32_t N = 0;

	uint64_t* minimizers = (uint64_t*)calloc(s_len, sizeof(uint64_t));
	for (int p = 0; p < w - 1; p++) {
		window[p] = hash32(GetSeedAtPos(s_2bit, p, k), k);
	}

	for (uint32_t p = w - 1; p < s_len - k - w; p++) {
		window[p%w] = hash32(GetSeedAtPos(s_2bit, p, k), k);
		m = Min_Window(window, w);
		if ((m != last_m) || (p - last_p >= w)) {
			minimizers[N++] = (m << 32) + p;
			last_m = m;
			last_p = p;
		}
	}

	uint64_t* new_minimizers = (uint64_t*)malloc(N * sizeof(uint64_t));
    std::memcpy(new_minimizers, minimizers, N * sizeof(uint64_t));
	delete[] minimizers;
	return std::pair <uint64_t*, uint32_t>(new_minimizers, N);
}

std::pair<uint64_t*, uint32_t> QTwoBitToMinimizers(uint32_t* s_2bit, uint32_t s_len, int k, int w) {
	uint32_t* window = (uint32_t*)calloc(w, sizeof(uint32_t));
	uint64_t last_m = 0;
	uint32_t last_p = 0;
	uint64_t m;

	uint32_t N = 0;

	uint64_t* minimizers = (uint64_t*)scalable_calloc(s_len, sizeof(uint64_t));
	for (int p = 0; p < w - 1; p++) {
		window[p] = hash32(GetSeedAtPos(s_2bit, p, k), k);
	}

	for (uint32_t p = w - 1; p < s_len - k; p++) {
		window[p%w] = hash32(GetSeedAtPos(s_2bit, p, k), k);
		m = Min_Window(window, w);
		if ((m != last_m) || (p - last_p >= w)) {
			minimizers[N++] = ((uint64_t)p << 32) + m;
			last_m = m;
			last_p = p;
		}
	}

	uint64_t* new_minimizers = (uint64_t*)scalable_malloc(N * sizeof(uint64_t));
    std::memcpy(new_minimizers, minimizers, N * sizeof(uint64_t));
	scalable_free(minimizers);
	return std::pair <uint64_t*, uint32_t>(new_minimizers, N);
}

