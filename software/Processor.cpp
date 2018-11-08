#include "Processor.h"
#include "tbb\scalable_allocator.h"

const char A_NT = 0;
const char C_NT = 1;
const char G_NT = 2;
const char T_NT = 3;
const char N_NT = 4;

static int sub_mat[25] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
static int gap_open = 0;
static int gap_extend = 0;
static int long_gap_open = 0;
static int long_gap_extend = 0;

int Nt2Int(char nt, int complement)
{
	int ret = N_NT;

	switch (nt) {
	case 'a':
	case 'A': ret = (complement) ? T_NT : A_NT;
		break;
	case 'c':
	case 'C': ret = (complement) ? G_NT : C_NT;
		break;
	case 'g':
	case 'G': ret = (complement) ? C_NT : G_NT;
		break;
	case 't':
	case 'T': ret = (complement) ? A_NT : T_NT;
		break;
	case 'n':
	case 'N': ret = N_NT;
		break;
	default:
		break;
	}

	return ret;
}

void InitializeScoringParams(size_t token, Darwin::AlignmentScoringParams& request, Darwin::AlignmentScoringParamsResponse& response)
{
	sub_mat[0] = request.sub_AA;
	sub_mat[1] = request.sub_AC;
	sub_mat[2] = request.sub_AG;
	sub_mat[3] = request.sub_AT;
	sub_mat[4] = request.sub_N;
	sub_mat[5] = request.sub_AC;
	sub_mat[6] = request.sub_CC;
	sub_mat[7] = request.sub_CG;
	sub_mat[8] = request.sub_CT;
	sub_mat[9] = request.sub_N;
	sub_mat[10] = request.sub_AG;
	sub_mat[11] = request.sub_CG;
	sub_mat[12] = request.sub_GG;
	sub_mat[13] = request.sub_GT;
	sub_mat[14] = request.sub_N;
	sub_mat[15] = request.sub_AT;
	sub_mat[16] = request.sub_CT;
	sub_mat[17] = request.sub_GT;
	sub_mat[18] = request.sub_TT;
	sub_mat[19] = request.sub_N;
	sub_mat[20] = request.sub_N;
	sub_mat[21] = request.sub_N;
	sub_mat[22] = request.sub_N;
	sub_mat[23] = request.sub_N;
	sub_mat[24] = request.sub_N;
	gap_open = request.gap_open;
	gap_extend = request.gap_extend;
	long_gap_open = request.long_gap_open;
	long_gap_extend = request.long_gap_extend;
	response.status = Darwin::Status::OK;
}

void InitializeMemory(size_t token, char* dram, Darwin::InitializeDRAMMessage& request, Darwin::InitializeDRAMMessageResponse& response)
{
}

void BatchAlignment(size_t token, char* dram, Darwin::BatchAlignmentInputFieldsDRAM& request, Darwin::BatchAlignmentResultDRAM& result)
{
	size_t num_requests = request.requests.size();
	int do_traceback = request.do_traceback;
	result.results.resize(num_requests);
	for (size_t r = 0; r < num_requests; r++) {
		Darwin::AlignmentInputFieldsDRAM req = request.requests[r];
		uint64_t ref_bases_start_addr = req.ref_bases_start_addr;
		uint64_t query_bases_start_addr = req.query_bases_start_addr;
		size_t ref_size = req.ref_size;
		size_t query_size = req.query_size;
		size_t max_tb_steps = req.max_tb_steps;

		uint8_t  align_fields = req.align_fields;
		int      reverse_ref = ((align_fields & (1 << 4)) >> 4);
		int      complement_ref = ((align_fields & (1 << 3)) >> 3);
		int      reverse_query = ((align_fields & (1 << 2)) >> 2);
		int      complement_query = ((align_fields & (1 << 1)) >> 1);
		int      start_end = (align_fields & 1);

		//assert(ref_size < MAX_TILE_SIZE);
		//assert(query_size < MAX_TILE_SIZE);

		size_t aligned_query_size = query_size + 1;
		size_t extra_query = aligned_query_size % 8;
		if (extra_query != 0)
			extra_query = 8 - extra_query;
		aligned_query_size += extra_query;

		size_t aligned_ref_size = ref_size + 1;
		size_t extra_ref = aligned_ref_size % 8;
		if (extra_ref != 0)
			extra_ref = 8 - extra_ref;
		aligned_ref_size += extra_ref;

		char* pMem = (char*)scalable_aligned_malloc(
			sizeof(int16_t) * 12 * aligned_query_size +
			sizeof(int16_t*) * aligned_ref_size +
			sizeof(int16_t) * aligned_ref_size * aligned_query_size,
			64);

		char* pM = pMem;

		int16_t* h_matrix_wr = (int16_t*)pM; pM += sizeof(int16_t) * aligned_query_size;
		int16_t* m_matrix_wr = (int16_t*)pM; pM += sizeof(int16_t) * aligned_query_size;
		int16_t* i_matrix_wr = (int16_t*)pM; pM += sizeof(int16_t) * aligned_query_size;
		int16_t* d_matrix_wr = (int16_t*)pM; pM += sizeof(int16_t) * aligned_query_size;
		int16_t* I_matrix_wr = (int16_t*)pM; pM += sizeof(int16_t) * aligned_query_size;
		int16_t* D_matrix_wr = (int16_t*)pM; pM += sizeof(int16_t) * aligned_query_size;

		int16_t* h_matrix_rd = (int16_t*)pM; pM += sizeof(int16_t) * aligned_query_size;
		int16_t* m_matrix_rd = (int16_t*)pM; pM += sizeof(int16_t) * aligned_query_size;
		int16_t* i_matrix_rd = (int16_t*)pM; pM += sizeof(int16_t) * aligned_query_size;
		int16_t* d_matrix_rd = (int16_t*)pM; pM += sizeof(int16_t) * aligned_query_size;
		int16_t* I_matrix_rd = (int16_t*)pM; pM += sizeof(int16_t) * aligned_query_size;
		int16_t* D_matrix_rd = (int16_t*)pM; pM += sizeof(int16_t) * aligned_query_size;

		int16_t** dir_matrix = (int16_t**)pM; pM += sizeof(int16_t*) * aligned_ref_size;

		for (int i = 0; i < ref_size + 1; i++) {
			dir_matrix[i] = (int16_t*)pM; pM += sizeof(int16_t) * aligned_query_size;
		}

		for (int i = 0; i < query_size + 1; i++) {
			h_matrix_rd[i] = 0;
			m_matrix_rd[i] = 0;
			i_matrix_rd[i] = -INF;
			d_matrix_rd[i] = -INF;
			I_matrix_rd[i] = -INF;
			D_matrix_rd[i] = -INF;

			h_matrix_wr[i] = 0;
			m_matrix_wr[i] = 0;
			i_matrix_wr[i] = -INF;
			d_matrix_wr[i] = -INF;
			I_matrix_wr[i] = -INF;
			D_matrix_wr[i] = -INF;
		}

		for (int i = 0; i < ref_size + 1; i++) {
			dir_matrix[i][0] = 0;
		}

		for (int j = 0; j < query_size + 1; j++) {
			dir_matrix[0][j] = 0;
		}

		int max_score = 0;
		int pos_score = 0;
		int max_i = 0;
		int max_j = 0;

		//const __m128i _mm_nt2Int_mask = _mm_set1_epi16(0x0f);

		for (int i = 1; i < ref_size + 1; i++) {
			for (int k = 1; k < aligned_query_size; k++) {
				m_matrix_rd[k] = m_matrix_wr[k];
				h_matrix_rd[k] = h_matrix_wr[k];
				i_matrix_rd[k] = i_matrix_wr[k];
				d_matrix_rd[k] = d_matrix_wr[k];
				I_matrix_rd[k] = I_matrix_wr[k];
				D_matrix_rd[k] = D_matrix_wr[k];
			}

			//j - row number; i - column number
			for (int j = 1; j < query_size + 1; j++) {
				int ref_nt = (reverse_ref) ? Nt2Int(dram[ref_bases_start_addr + ref_size - i], complement_ref) : Nt2Int(dram[ref_bases_start_addr + i - 1], complement_ref);
				int query_nt = (reverse_query) ? Nt2Int(dram[query_bases_start_addr + query_size - j], complement_query) : Nt2Int(dram[query_bases_start_addr + j - 1], complement_query);

				int match;

				//value from the W matrix for the match/mismatch penalty/point
				match = sub_mat[query_nt * 5 + ref_nt];

                m_matrix_wr[j] = h_matrix_rd[j - 1] + match; 

				int del_open = m_matrix_rd[j] + gap_open;
				int del_extend = d_matrix_rd[j] + gap_extend;
				int ins_open = m_matrix_wr[j - 1] + gap_open;
				int ins_extend = i_matrix_wr[j - 1] + gap_extend;

				int long_del_open = m_matrix_rd[j] + long_gap_open;
				int long_del_extend = D_matrix_rd[j] + long_gap_extend;
				int long_ins_open = m_matrix_wr[j - 1] + long_gap_open;
				int long_ins_extend = I_matrix_wr[j - 1] + long_gap_extend;

				(dir_matrix)[i][j] = 0;
                if (ins_open > ins_extend) {
                    i_matrix_wr[j] = ins_open;
                    (dir_matrix)[i][j] += (8 << INSERT_OP);
                }
                else {
                    i_matrix_wr[j] = ins_extend;
                    (dir_matrix)[i][j] += 0;
                }

                if (del_open > del_extend) {
                    d_matrix_wr[j] = del_open;
                    (dir_matrix)[i][j] += (8 << DELETE_OP);
                }
                else {
                    d_matrix_wr[j] = del_extend;
                    (dir_matrix)[i][j] += 0;
                }
                
                if (long_ins_open > long_ins_extend) {
                    I_matrix_wr[j] = long_ins_open;
                    (dir_matrix)[i][j] += (8 << LONG_INSERT_OP);
                }
                else {
                    I_matrix_wr[j] = long_ins_extend;
                    (dir_matrix)[i][j] += 0;
                }

                if (long_del_open > long_del_extend) {
                    D_matrix_wr[j] = long_del_open;
                    (dir_matrix)[i][j] += (8 << LONG_DELETE_OP);
                }
                else {
                    D_matrix_wr[j] = long_del_extend;
                    (dir_matrix)[i][j] += 0;
                }


                if ((m_matrix_wr[j] >= 0) && (m_matrix_wr[j] >= d_matrix_wr[j]) && (m_matrix_wr[j] >= D_matrix_wr[j]) && (m_matrix_wr[j] >= i_matrix_wr[j]) && (m_matrix_wr[j] >= I_matrix_wr[j])) {
                    h_matrix_wr[j] = m_matrix_wr[j];
                    (dir_matrix)[i][j] += M;
                }
                else if ((d_matrix_wr[j] >= 0) && (d_matrix_wr[j] >= D_matrix_wr[j]) && (d_matrix_wr[j] >= i_matrix_wr[j]) && (d_matrix_wr[j] >= I_matrix_wr[j])) {
                    h_matrix_wr[j] = d_matrix_wr[j];
                    (dir_matrix)[i][j] += D;
                }
                else if ((D_matrix_wr[j] >= 0) && (D_matrix_wr[j] >= i_matrix_wr[j]) && (D_matrix_wr[j] >= I_matrix_wr[j])) {
                    h_matrix_wr[j] = D_matrix_wr[j];
                    (dir_matrix)[i][j] += L_D;
                }
                else if ((i_matrix_wr[j] >= 0) && (i_matrix_wr[j] >= I_matrix_wr[j])) {
                    h_matrix_wr[j] = i_matrix_wr[j];
                    (dir_matrix)[i][j] += I;
                }
                else if ((I_matrix_wr[j] >= 0)) { 
                    h_matrix_wr[j] = I_matrix_wr[j];
                    (dir_matrix)[i][j] += L_I;
                }
                else {
                    h_matrix_wr[j] = 0;
                    (dir_matrix)[i][j] += Z;
                }

				//std::cerr << i_matrix_wr[j] << " " << d_matrix_wr[j] << " " << m_matrix_wr[j] << " " << I_matrix_wr[j] << " " << D_matrix_wr[j] << " " << ((dir_matrix)[i][j]) << std::endl;

				if (h_matrix_wr[j] >= max_score) {
					max_score = h_matrix_wr[j];
					max_i = i;
					max_j = j;
				}

				pos_score = h_matrix_wr[j];

			}
		}


		result.results[r].index = req.index;
		if (start_end) {
			result.results[r].ref_max_pos = ref_size - 1;
			result.results[r].query_max_pos = query_size - 1;
			result.results[r].score = pos_score;
		}
		else {
			result.results[r].ref_max_pos = max_i - 1;
			result.results[r].query_max_pos = max_j - 1;
			result.results[r].score = max_score;
		}
		uint64_t tb_pointer = 0;
		int      total_TB_pointers = 0;
		result.results[r].TB_pointers.clear();
		if (do_traceback) {
			int i_curr = max_i, j_curr = max_j;
			if (start_end) {
				i_curr = ref_size;
				j_curr = query_size;
			}
			int i_steps = 0, j_steps = 0;
			uint64_t state = dir_matrix[i_curr][j_curr] % 8;

			while (state != Z) {
				if ((i_curr < 0) || (j_curr < 0)) {
					break;
				}
				if ((i_steps == max_tb_steps) || (j_steps == max_tb_steps)) {
					break;
				}
                
                uint64_t tb_state = state % 4;
				if (total_TB_pointers % 32 == 0) {
					if (total_TB_pointers > 0) {
						result.results[r].TB_pointers.push_back(tb_pointer);
					}
					tb_pointer = tb_state;
				}
				else {
					tb_pointer = (tb_state << (2 * (total_TB_pointers % 32))) + tb_pointer;
				}
				total_TB_pointers++;

				if (state == M) {
                    if ((i_curr > 0) && (j_curr > 0))
                        state = (dir_matrix[i_curr - 1][j_curr - 1] % 8);
                    else
                        state = Z;
					i_curr--;
					j_curr--;
					i_steps++;
					j_steps++;
				}
				else if (state == I) {
					state = (dir_matrix[i_curr][j_curr] & (8 << INSERT_OP)) ? M : I;
					j_curr--;
					j_steps++;
				}
				else if (state == D) {
					state = (dir_matrix[i_curr][j_curr] & (8 << DELETE_OP)) ? M : D;
					i_curr--;
					i_steps++;
				}
				else if (state == L_I) {
					state = (dir_matrix[i_curr][j_curr] & (8 << LONG_INSERT_OP)) ? M : L_I;
					j_curr--;
					j_steps++;
				}
				else if (state == L_D) {
					state = (dir_matrix[i_curr][j_curr] & (8 << LONG_DELETE_OP)) ? M : L_D;
					i_curr--;
					i_steps++;
				}

			}

			if (total_TB_pointers > 0) {
				result.results[r].TB_pointers.push_back(tb_pointer);
			}

			result.results[r].ref_offset = i_steps;
			result.results[r].query_offset = j_steps;
			result.results[r].total_TB_pointers = total_TB_pointers;

		}
		else {
			result.results[r].ref_offset = 0;
			result.results[r].query_offset = 0;
			result.results[r].total_TB_pointers = 0;
		}

		scalable_aligned_free(pMem);
	}
	std::random_shuffle(result.results.begin(), result.results.end());
}

DRAM *g_DRAM = nullptr;

InitializeScoringParameters_ptr g_InitializeScoringParameters = InitializeScoringParams;
InitializeReferenceMemory_ptr g_InitializeReferenceMemory = InitializeMemory; 
InitializeReadMemory_ptr g_InitializeReadMemory = InitializeMemory;

BatchAlignment_ptr g_BatchAlignment = BatchAlignment;
//BatchAlignmentTwoPiece_ptr g_BatchAlignmentTwoPiece = BatchAlignmentTwoPiece;
