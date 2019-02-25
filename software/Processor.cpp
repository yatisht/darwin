#include "Processor.h"

#define SWAP(A,B) { __m256i* tmp = A; A = B; B = tmp; }
#define SWAP3(A,B,C) { __m256i* tmp = A; A = B; B = C; C = tmp; }
#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)

const char A_NT = 0;
const char C_NT = 1;
const char G_NT = 2;
const char T_NT = 3;
const char N_NT = 4;
const int NUM_AMINO_ACIDS = 5; // including N
const int16_t NEG_INF = (INT16_MIN/(int16_t)(2));

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

__m256i* CreateVProfile(
		const char * const dram, uint64_t query_bases_start_addr, const int queryLen, bool reverse_query, bool complement)
{
	int32_t i = 0;
	int32_t j = 0;
	int32_t k = 0;
	int32_t segNum = 0;
	const int32_t n = NUM_AMINO_ACIDS; /* number of amino acids in table */
	const int32_t segWidth = 16; /* number of values in vector unit */
	const int32_t segLen = (queryLen + segWidth - 1) / segWidth;
	__m256i* vProfile = (__m256i *) scalable_aligned_malloc(sizeof(__m256i) * n * segLen, 32);
	int32_t index = 0;

	for (k=0; k<n; ++k) {
		for (i=0; i<segLen; ++i) {
			__m256i_16_t t;
			j = i;
			for (segNum=0; segNum<segWidth; ++segNum) {
				uint64_t char_idx = reverse_query ? query_bases_start_addr + (queryLen - 1) - j : query_bases_start_addr + j;
				t.v[segNum] = j >= queryLen ? 0 : sub_mat[n*k+Nt2Int(dram[char_idx], complement)];
				j += segLen;
			}
			_mm256_store_si256(&vProfile[index], t.m);
			++index;
		}
	}

	return vProfile;
}

result_t* result_new_trace(const int a, const int b, const size_t alignment, const size_t size)
{
	/* declare all variables */
	result_t *result = (result_t*)scalable_malloc(sizeof(result_t));
	result->score = 0;
	result->end_query = 0;
	result->end_ref = 0;
	result->trace = (__m256i *) scalable_aligned_malloc(size*a*b*sizeof(__m256i), alignment);
	return result;
}

void memset___m256i(__m256i *b, __m256i c, size_t len)
{
	size_t i;
	for (i=0; i<len; ++i) {
		_mm256_store_si256(&b[i], c);
	}
}

static inline int16_t _mm256_hmax_epi16_rpl(__m256i a) {
	a = _mm256_max_epi16(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,0,0)));
	a = _mm256_max_epi16(a, _mm256_slli_si256(a, 8));
	a = _mm256_max_epi16(a, _mm256_slli_si256(a, 4));
	a = _mm256_max_epi16(a, _mm256_slli_si256(a, 2));
	return _mm256_extract_epi16(a, 15);
}


static inline void arr_store(
		__m256i *array,
		__m256i vH,
		int32_t t,
		int32_t seglen,
		int32_t d)
{
	_mm256_store_si256(array + (1LL*d*seglen+t), vH);
}

static inline __m256i arr_load(
		__m256i *array,
		int32_t t,
		int32_t seglen,
		int32_t d)
{
	return _mm256_load_si256(array + (1LL*d*seglen+t));
}

result_t* DualAlignSIMD(
		const char * const dram,
		uint64_t query_bases_start_addr, const int queryLen,
		uint64_t ref_bases_start_addr, const int refLen,
		bool reverse_query, bool reverse_ref,
		bool complement_query, bool complement_ref,
		bool start_end)
{

	const int32_t segWidth = 16; /* number of values in vector unit */
	const int32_t segLen = (queryLen + segWidth - 1) / segWidth;

	result_t *result = result_new_trace(segLen, refLen, 32, sizeof(__m256i));
	if (queryLen == 0 || refLen == 0) {
		result->score = 0;
		result->end_query = 0;
		result->end_ref = 0;
		return result;
	}

    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;

    __m256i* const vProfile = CreateVProfile(dram, query_bases_start_addr, queryLen, reverse_query, complement_query);
    __m256i* pvHStore = (__m256i *) scalable_aligned_malloc(segLen * sizeof(__m256i), 32);
    __m256i* pvHLoad = (__m256i *) scalable_aligned_malloc(segLen * sizeof(__m256i), 32);
    __m256i* const pvE = (__m256i *) scalable_aligned_malloc(segLen * sizeof(__m256i), 32);
    __m256i* pvEaStore = (__m256i *) scalable_aligned_malloc(segLen * sizeof(__m256i), 32);
    __m256i* pvEaLoad = (__m256i *) scalable_aligned_malloc(segLen * sizeof(__m256i), 32);
    __m256i* const pvE_L = (__m256i *) scalable_aligned_malloc(sizeof(__m256i) * segLen, 32);
    __m256i* pvE_LaStore = (__m256i *) scalable_aligned_malloc(sizeof(__m256i) * segLen, 32);
    __m256i* pvE_LaLoad = (__m256i *) scalable_aligned_malloc(sizeof(__m256i) * segLen, 32);
    
    __m256i* const pvHT = (__m256i *) scalable_aligned_malloc(segLen * sizeof(__m256i), 32);
    __m256i* pvHMax = (__m256i *) scalable_aligned_malloc(segLen * sizeof(__m256i), 32);
    __m256i vGapO = _mm256_set1_epi16(-gap_open);
    __m256i vGapE = _mm256_set1_epi16(-gap_extend);
    __m256i vGapO_L = _mm256_set1_epi16(-long_gap_open);
    __m256i vGapE_L = _mm256_set1_epi16(-long_gap_extend);
    __m256i vZero = _mm256_setzero_si256();
    int16_t score = 0;
    __m256i vMaxH = vZero;
    __m256i vMaxHUnit = vZero;
    __m256i vTZeroT = _mm256_set1_epi16(TRACEBACK_ZERO_T);
    __m256i vTDelT  = _mm256_set1_epi16(TRACEBACK_DEL_T);
    __m256i vTInsT  = _mm256_set1_epi16(TRACEBACK_INS_T);
    __m256i vTDel_LT = _mm256_set1_epi16(TRACEBACK_DEL_L_T);
    __m256i vTIns_LT = _mm256_set1_epi16(TRACEBACK_INS_L_T);
    __m256i vTDiag = _mm256_set1_epi16(TRACEBACK_DIAG);
    __m256i vTDiagDel = _mm256_set1_epi16(TRACEBACK_DIAG_DEL);
    __m256i vTDel = _mm256_set1_epi16(TRACEBACK_DEL);
    __m256i vTDiagIns = _mm256_set1_epi16(TRACEBACK_DIAG_INS);
    __m256i vTIns = _mm256_set1_epi16(TRACEBACK_INS);
    __m256i vTDiagDel_L = _mm256_set1_epi16(TRACEBACK_DIAG_DEL_L);
    __m256i vTDel_L = _mm256_set1_epi16(TRACEBACK_DEL_L);
    __m256i vTDiagIns_L = _mm256_set1_epi16(TRACEBACK_DIAG_INS_L);
    __m256i vTIns_L = _mm256_set1_epi16(TRACEBACK_INS_L);
    __m256i vTMask = _mm256_set1_epi16(TRACEBACK_T_MASK);
    __m256i vFTMask = _mm256_set1_epi16(TRACEBACK_F_MASK);
    __m256i vF_LTMask = _mm256_set1_epi16(TRACEBACK_F_L_MASK);

    /* initialize H and E */
    memset___m256i(pvHStore, vZero, segLen);
    memset___m256i(pvE, _mm256_set1_epi16(gap_open), segLen);
    memset___m256i(pvEaStore, _mm256_set1_epi16(gap_open), segLen);

	// Initialize E_L
	memset___m256i(pvE_L, _mm256_set1_epi16(long_gap_open), segLen);
	memset___m256i(pvE_LaStore, _mm256_set1_epi16(long_gap_open), segLen);

	for (i=0; i<segLen; ++i) {
		arr_store(result->trace, _mm256_add_epi16(vTDiagDel, vTDiagDel_L), i, segLen, 0);
	}

	/* outer loop over database sequence */
	for (j=0; j<refLen; ++j) {
		__m256i vEF_opn;
		__m256i vE;
		__m256i vE_ext;
		__m256i vF;
		__m256i vF_ext;
		__m256i vFa;
		__m256i vFa_ext;
		__m256i vH;
		__m256i vH_dag;

		__m256i vEF_L_opn;
		__m256i vE_L;
		__m256i vE_L_ext;

		__m256i vF_L;
		__m256i vF_L_ext;
		__m256i vF_La;
		__m256i vF_La_ext;
		const __m256i* vP = NULL;

        /* Initialize F value to 0.  Any errors to vH values will be
         * corrected in the Lazy_F loop. */
        vF = _mm256_sub_epi16(vZero,vGapO);

		/* Initialize F value to 0.  Any errors to vH values will be
         * corrected in the Lazy_F loop. */
		vF_L = _mm256_sub_epi16(vZero,vGapO_L);

        /* load final segment of pvHStore and shift left by 2 bytes */
        vH = _mm256_load_si256(&pvHStore[segLen - 1]);
        vH = _mm256_slli_si256_rpl(vH, 2);

		/* Correct part of the vProfile */
		uint64_t char_idx = reverse_ref ? ref_bases_start_addr + (refLen - 1) - j : ref_bases_start_addr + j;
		vP = vProfile + Nt2Int(dram[char_idx], complement_ref) * segLen; // TODO: complement

		if (end_ref == j-2) {
			/* Swap in the max buffer. */
			SWAP3(pvHMax,  pvHLoad,  pvHStore)
			SWAP(pvEaLoad,  pvEaStore)
			SWAP(pvE_LaLoad, pvE_LaStore)
		}
		else {
			/* Swap the 2 H buffers. */
			SWAP(pvHLoad,  pvHStore)
			SWAP(pvEaLoad,  pvEaStore)
			SWAP(pvE_LaLoad, pvE_LaStore)
		}

		/* inner loop to process the query sequence */
		for (i=0; i<segLen; ++i) {
			vE = _mm256_load_si256(pvE + i);
			vE_L = _mm256_load_si256(pvE_L + i);

			/* Get max from vH, vE and vF. */
			vH_dag = _mm256_add_epi16(vH, _mm256_load_si256(vP + i));
			vH_dag = _mm256_max_epi16(vH_dag, vZero);
			vH = _mm256_max_epi16(vH_dag, vE);
			vH = _mm256_max_epi16(vH, vF);
			vH = _mm256_max_epi16(vH, vE_L);
			vH = _mm256_max_epi16(vH, vF_L);
			/* Save vH values. */
			_mm256_store_si256(pvHStore + i, vH);

			{
				__m256i vTAll = arr_load(result->trace, i, segLen, j);
				__m256i cond_zero = _mm256_cmpeq_epi16(vH, vZero); // is zero
				__m256i case1 = _mm256_cmpeq_epi16(vH, vH_dag);    // is match
				__m256i case2 = _mm256_cmpeq_epi16(vH, vF);        // is delete
				__m256i case3 = _mm256_cmpeq_epi16(vH, vF_L);      // is delete_long
				__m256i case4 = _mm256_cmpeq_epi16(vH, vE_L);      // is insert_long
				__m256i vT = _mm256_blendv_epi8(
						_mm256_blendv_epi8(vTDelT, vTInsT, case2),
						_mm256_blendv_epi8(vTDiag, vTZeroT, cond_zero),
						case1);
				vT = _mm256_blendv_epi8(
						vT,
						vTIns_LT,
						case3);
				vT = _mm256_blendv_epi8(
						vT,
						vTDel_LT,
						case4);
				_mm256_store_si256(pvHT + i, vT);
				vT = _mm256_or_si256(vT, vTAll);
				arr_store(result->trace, vT, i, segLen, j);
			}

			vMaxH = _mm256_max_epi16(vH, vMaxH);
			vEF_opn = _mm256_sub_epi16(vH, vGapO);
			vEF_L_opn = _mm256_sub_epi16(vH, vGapO_L);

			/* Update vE/E_L value. */
			vE_ext = _mm256_sub_epi16(vE, vGapE);
			vE = _mm256_max_epi16(vEF_opn, vE_ext);
			_mm256_store_si256(pvE + i, vE);
			vE_L_ext = _mm256_sub_epi16(vE_L, vGapE_L);
			vE_L = _mm256_max_epi16(vEF_L_opn, vE_L_ext);
			_mm256_store_si256(pvE_L + i, vE_L);
			{
				__m256i vEa = _mm256_load_si256(pvEaLoad + i);
				__m256i vEa_ext = _mm256_sub_epi16(vEa, vGapE);
				__m256i vE_La = _mm256_load_si256(pvE_LaLoad + i);
				__m256i vE_La_ext = _mm256_sub_epi16(vE_La, vGapE_L);
				vEa = _mm256_max_epi16(vEF_opn, vEa_ext);
				_mm256_store_si256(pvEaStore + i, vEa);
				vE_La = _mm256_max_epi16(vEF_L_opn, vE_La_ext);
				_mm256_store_si256(pvE_LaStore + i, vE_La);

				if (j+1<refLen) {
					__m256i cond_E = _mm256_cmpgt_epi16(vEF_opn, vEa_ext);
					__m256i vT_E = _mm256_blendv_epi8(vTDel, vTDiagDel, cond_E);
					__m256i cond_E_L = _mm256_cmpgt_epi16(vEF_L_opn, vE_La_ext);
					__m256i vT_E_L = _mm256_blendv_epi8(vTDel_L, vTDiagDel_L, cond_E_L);
					__m256i vT = _mm256_or_si256(vT_E, vT_E_L);
					arr_store(result->trace, vT, i, segLen, j+1);
				}
			}

			/* Update vF/vF_L value. */
			vF_ext = _mm256_sub_epi16(vF, vGapE);
			vF = _mm256_max_epi16(vEF_opn, vF_ext);
			vF_L_ext = _mm256_sub_epi16(vF_L, vGapE_L);
			vF_L = _mm256_max_epi16(vEF_L_opn, vF_L_ext);
			if (i+1<segLen) {
				__m256i vTAll = arr_load(result->trace, i+1, segLen, j);
				__m256i cond_F = _mm256_cmpgt_epi16(vEF_opn, vF_ext);
				__m256i vT_F = _mm256_blendv_epi8(vTIns, vTDiagIns, cond_F);
				__m256i cond_F_L = _mm256_cmpgt_epi16(vEF_L_opn, vF_L_ext);
				__m256i vT_F_L = _mm256_blendv_epi8(vTIns_L, vTDiagIns_L, cond_F_L);
				__m256i vT = _mm256_or_si256(vT_F, vT_F_L);
				vT = _mm256_or_si256(vT, vTAll);
				arr_store(result->trace, vT, i+1, segLen, j);
			}

			/* Load the next vH. */
			vH = _mm256_load_si256(pvHLoad + i);
		}

		/* Lazy_F/F_L loop */
        vFa_ext = vF_ext;
        vFa = vF;
        for (k=0; k<segWidth; ++k) {
            __m256i vHp = _mm256_load_si256(&pvHLoad[segLen - 1]); // load the last value added to pvHLoad
            vHp = _mm256_slli_si256_rpl(vHp, 2);
            vEF_opn = _mm256_slli_si256_rpl(vEF_opn, 2);
            vEF_opn = _mm256_insert_epi16(vEF_opn, gap_open, 0);
            vF_ext = _mm256_slli_si256_rpl(vF_ext, 2);
            vF_ext = _mm256_insert_epi16(vF_ext, NEG_INF, 0);
            vF = _mm256_slli_si256_rpl(vF, 2);
            vF = _mm256_insert_epi16(vF, gap_open, 0);
            vFa_ext = _mm256_slli_si256_rpl(vFa_ext, 2);
            vFa_ext = _mm256_insert_epi16(vFa_ext, NEG_INF, 0);
            vFa = _mm256_slli_si256_rpl(vFa, 2);
            vFa = _mm256_insert_epi16(vFa, gap_open, 0);

			vEF_L_opn = _mm256_slli_si256_rpl(vEF_L_opn, 2);
			vEF_L_opn = _mm256_insert_epi16(vEF_L_opn, long_gap_open, 0);
			vF_L_ext = _mm256_slli_si256_rpl(vF_L_ext, 2);
			vF_L_ext = _mm256_insert_epi16(vF_L_ext, NEG_INF, 0);
			vF_L = _mm256_slli_si256_rpl(vF_L, 2);
			vF_L = _mm256_insert_epi16(vF_L, long_gap_open, 0);
			vF_La_ext = _mm256_slli_si256_rpl(vF_La_ext, 2);
			vF_La_ext = _mm256_insert_epi16(vF_La_ext, NEG_INF, 0);
			vF_La = _mm256_slli_si256_rpl(vF_La, 2);
			vF_La = _mm256_insert_epi16(vF_La, long_gap_open, 0);
            for (i=0; i<segLen; ++i) {
                vH = _mm256_load_si256(pvHStore + i);
                vH = _mm256_max_epi16(vH,vF);
                vH = _mm256_max_epi16(vH,vF_L);
                _mm256_store_si256(pvHStore + i, vH);
                {
                    __m256i vTAll;
                    __m256i vT;
                    __m256i case1;
                    __m256i case2;
                    __m256i case3;
                    __m256i cond_F;
                    __m256i cond_F_L;
                    vHp = _mm256_add_epi16(vHp, _mm256_load_si256(vP + i));
                    vHp = _mm256_max_epi16(vHp, vZero);
                    case1 = _mm256_cmpeq_epi16(vH, vHp);          // vH = vHp
                    case2 = _mm256_cmpeq_epi16(vH, vF);           // vH = vF
                    case3 = _mm256_cmpeq_epi16(vH, vF_L);         // vH = vF_L
                    cond_F = _mm256_andnot_si256(case1,case2);    // !vHp & vF
                    cond_F_L = _mm256_andnot_si256(_mm256_or_si256(case1, case2), case3); // !(vHp | vF) & vF_L
                    vT = _mm256_load_si256(pvHT + i);
                    vT = _mm256_blendv_epi8(vT, vTInsT, cond_F);
					vT = _mm256_blendv_epi8(vT, vTIns_LT, cond_F_L);
                    _mm256_store_si256(pvHT + i, vT);
                    vTAll = arr_load(result->trace, i, segLen, j);
                    vTAll = _mm256_and_si256(vTAll, vTMask);
                    vTAll = _mm256_or_si256(vTAll, vT);
                    arr_store(result->trace, vTAll, i, segLen, j);
                }
                vMaxH = _mm256_max_epi16(vH, vMaxH);
                /* Update vF value. */
                {
                    __m256i vTAll = arr_load(result->trace, i, segLen, j);
                    __m256i cond_F = _mm256_cmpgt_epi16(vEF_opn, vFa_ext);
                    __m256i vT_F = _mm256_blendv_epi8(vTIns, vTDiagIns, cond_F);
					__m256i cond_F_L = _mm256_cmpgt_epi16(vEF_L_opn, vF_La_ext);
					__m256i vT_F_L = _mm256_blendv_epi8(vTIns_L, vTDiagIns_L, cond_F_L);
                    vTAll = _mm256_and_si256(vTAll, vFTMask);
                    vTAll = _mm256_or_si256(vTAll, vT_F);
					vTAll = _mm256_and_si256(vTAll, vF_LTMask);
					vTAll = _mm256_or_si256(vTAll, vT_F_L);
                    arr_store(result->trace, vTAll, i, segLen, j);
                }
				/* Update vE/E_L value. */
                vEF_opn = _mm256_sub_epi16(vH, vGapO);
                vF_ext = _mm256_sub_epi16(vF, vGapE);
				vEF_L_opn = _mm256_sub_epi16(vH, vGapO_L);
				vF_L_ext = _mm256_sub_epi16(vF_L, vGapE_L);
                {
                    __m256i vEa = _mm256_load_si256(pvEaLoad + i);
                    __m256i vEa_ext = _mm256_sub_epi16(vEa, vGapE);
                    vEa = _mm256_max_epi16(vEF_opn, vEa_ext);
                    _mm256_store_si256(pvEaStore + i, vEa);

					__m256i vE_La = _mm256_load_si256(pvE_LaLoad + i);
					__m256i vE_La_ext = _mm256_sub_epi16(vE_La, vGapE_L);
					vE_La = _mm256_max_epi16(vEF_L_opn, vE_La_ext);
					_mm256_store_si256(pvE_LaStore + i, vE_La);
                    if (j+1<refLen) {
                        __m256i cond_E = _mm256_cmpgt_epi16(vEF_opn, vEa_ext);
                        __m256i vT_E = _mm256_blendv_epi8(vTDel, vTDiagDel, cond_E);
						__m256i cond_E_L = _mm256_cmpgt_epi16(vEF_L_opn, vE_La_ext);
						__m256i vT_E_L = _mm256_blendv_epi8(vTDel_L, vTDiagDel_L, cond_E_L);
						__m256i vT = _mm256_or_si256(vT_E, vT_E_L);
                        arr_store(result->trace, vT, i, segLen, j+1);
                    }
                }
                if (! _mm256_movemask_epi8(
                        _mm256_or_si256(
                                _mm256_cmpgt_epi16(vF_ext, vEF_opn),
                                _mm256_cmpeq_epi16(vF_ext, vEF_opn))) &&
					! _mm256_movemask_epi8(
							_mm256_or_si256(
									_mm256_cmpgt_epi16(vF_L_ext, vEF_L_opn),
									_mm256_cmpeq_epi16(vF_L_ext, vEF_L_opn)))

                                )
                    goto end;
                /*vF = _mm256_max_epi16(vEF_opn, vF_ext);*/
                vF = vF_ext;
                vFa_ext = _mm256_sub_epi16(vFa, vGapE);
                vFa = _mm256_max_epi16(vEF_opn, vFa_ext);

				vF_L = vF_L_ext;
				vF_La_ext = _mm256_sub_epi16(vF_La, vGapE_L);
				vF_La = _mm256_max_epi16(vEF_L_opn, vF_La_ext);
				vHp = _mm256_load_si256(pvHLoad + i);
            }
        }
        end:
        {
        }

        {
            __m256i vCompare = _mm256_cmpgt_epi16(vMaxH, vMaxHUnit);
            if (_mm256_movemask_epi8(vCompare)) {
                score = _mm256_hmax_epi16_rpl(vMaxH);
                vMaxHUnit = _mm256_set1_epi16(score);
                end_ref = j;
            }
        }

        /*if (score == stop) break;*/
    }

	if (start_end) {
		/* if start_end extract final vH */
		__m256i_16_t* end_segment = (__m256i_16_t *) &pvHStore[(queryLen - 1) % segLen];
		score = end_segment->v[(queryLen - 1) / segLen];
	} else {
		/* otherwise extract max */
		if (end_ref == j-1) {
			/* end_ref was the last store column */
			SWAP(pvHMax,  pvHStore)
		}
		else if (end_ref == j-2) {
			/* end_ref was the last load column */
			SWAP(pvHMax,  pvHLoad)
		}
		/* Trace the alignment ending position on read. */
		{
			int16_t *t = (int16_t *) pvHMax;
			int32_t column_len = segLen * segWidth;
			end_query = queryLen - 1;
			for (i = 0; i < column_len; ++i, ++t) {
				if (*t == score) {
					int32_t temp = i / segWidth + i % segWidth * segLen;
					if (temp < end_query) {
						end_query = temp;
					}
				}
			}
		}
	}

	if (start_end) {
		result->score = score;
		result->end_query = queryLen - 1;
		result->end_ref = refLen - 1;
	} else {
		result->score = score;
		result->end_query = end_query;
		result->end_ref = end_ref;
	}

    scalable_aligned_free(pvHMax);
    scalable_aligned_free(pvHT);
	scalable_aligned_free(pvE_LaLoad);
	scalable_aligned_free(pvE_LaStore);
	scalable_aligned_free(pvE_L);
    scalable_aligned_free(pvEaLoad);
    scalable_aligned_free(pvEaStore);
    scalable_aligned_free(pvE);
    scalable_aligned_free(pvHLoad);
    scalable_aligned_free(pvHStore);

    return result;
}

void AddToTracebackPointers(uint64_t state, AlignmentResult &alignmentResult, uint64_t &tb_pointer, uint16_t &total_TB_pointers) {

	uint64_t tb_state = state % 4;
	if (total_TB_pointers % 32 == 0) {
		if (total_TB_pointers > 0) {
			alignmentResult.TB_pointers.push_back(tb_pointer);
		}
		tb_pointer = tb_state;
	}
	else {
		tb_pointer = (tb_state << (2 * (total_TB_pointers % 32))) + tb_pointer;
	}
	total_TB_pointers++;

}

// TODO: start end score
void DualAlignSIMDTraceback(
		result_t *result,
		int len_query,
		int len_ref,
		int max_tb_steps,
		AlignmentResult &alignmentResult,
		bool start_end) {

	int64_t i = result->end_query;
	int64_t j = result->end_ref;
	if (start_end) {
		i = len_query - 1;
		j = len_ref - 1;
	}

	assert(i < len_query);
	assert(j < len_ref);

	uint16_t i_steps = 0, j_steps = 0;
	int where = TRACEBACK_DIAG;
	int16_t *HT = (int16_t*)result->trace;
	int64_t segWidth = 16;
	int64_t segLen = 0;
	segLen = (len_query + segWidth - 1) / segWidth;

	uint64_t tb_pointer = 0;
	uint16_t total_TB_pointers = 0;
	alignmentResult.TB_pointers.clear();
	while (i >= 0 && j >= 0) {
		int64_t loc = j*segLen*segWidth + (i%segLen)*segWidth + (i/segLen);

		if ((i_steps == max_tb_steps) || (j_steps == max_tb_steps)) {
			break;
		}

		if (TRACEBACK_DIAG == where) {
			if (HT[loc] & TRACEBACK_DIAG) {
				AddToTracebackPointers(M, alignmentResult, tb_pointer, total_TB_pointers);
				--i;
				--j;
				i_steps++;
				j_steps++;
			}
			else if (HT[loc] & TRACEBACK_DEL_T) {
				where = TRACEBACK_DEL_T;
			}
			else if (HT[loc] & TRACEBACK_INS_T) {
				where = TRACEBACK_INS_T;
			}
			else if (HT[loc] & TRACEBACK_DEL_L_T) {
				where = TRACEBACK_DEL_L_T;
			}
			else if (HT[loc] & TRACEBACK_INS_L_T) {
				where = TRACEBACK_INS_L_T;
			}
			else {
				break;
			}
		}
		else if (TRACEBACK_DEL_T == where) {
			AddToTracebackPointers(D, alignmentResult, tb_pointer, total_TB_pointers);
			--j;
			j_steps++;
			if (HT[loc] & TRACEBACK_DIAG_DEL) {
				where = TRACEBACK_DIAG;
			}
			else if (HT[loc] & TRACEBACK_DEL) {
				where = TRACEBACK_DEL_T;
			}
			else {
				assert(0);
			}
		}
		else if (TRACEBACK_INS_T == where) {
			AddToTracebackPointers(I, alignmentResult, tb_pointer, total_TB_pointers);
			--i;
			i_steps++;
			if (HT[loc] & TRACEBACK_DIAG_INS) {
				where = TRACEBACK_DIAG;
			}
			else if (HT[loc] & TRACEBACK_INS) {
				where = TRACEBACK_INS_T;
			}
			else {
				assert(0);
			}
		}
		else if (TRACEBACK_INS_L_T == where) {
			AddToTracebackPointers(L_I, alignmentResult, tb_pointer, total_TB_pointers);
			--i;
			i_steps++;
			if (HT[loc] & TRACEBACK_DIAG_INS_L) {
				where = TRACEBACK_DIAG;
			}
			else if (HT[loc] & TRACEBACK_INS_L) {
				where = TRACEBACK_INS_L_T;
			}
			else {
				assert(0);
			}
		}
		else if (TRACEBACK_DEL_L_T == where) {
			AddToTracebackPointers(L_D, alignmentResult, tb_pointer, total_TB_pointers);
			--j;
			j_steps++;
			if (HT[loc] & TRACEBACK_DIAG_DEL_L) {
				where = TRACEBACK_DIAG;
			}
			else if (HT[loc] & TRACEBACK_DEL_L) {
				where = TRACEBACK_DEL_L_T;
			}
			else {
				assert(0);
			}

		}
		else if (TRACEBACK_ZERO_T == where) {
			break;
		}
		else {
			assert(0);
		}
	}

	if (total_TB_pointers > 0) {
		alignmentResult.TB_pointers.push_back(tb_pointer);
	}

	alignmentResult.query_offset = i_steps;
	alignmentResult.ref_offset = j_steps;
	alignmentResult.total_TB_pointers = total_TB_pointers;
}

void BatchAlignmentSIMD(size_t token, char* dram, Darwin::BatchAlignmentInputFieldsDRAM& request, Darwin::BatchAlignmentResultDRAM& result) {
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
		int 	 reverse_ref = ((align_fields & (1 << 4)) >> 4);
		int      complement_ref = ((align_fields & (1 << 3)) >> 3);
		int      reverse_query = ((align_fields & (1 << 2)) >> 2);
		int      complement_query = ((align_fields & (1 << 1)) >> 1);
		int      start_end = (align_fields & 1);

		result_t *simd_result = DualAlignSIMD(dram,
				query_bases_start_addr, query_size,
				ref_bases_start_addr, ref_size,
				reverse_query, reverse_ref,
				complement_query, complement_ref,
				start_end);

		result.results[r].index = req.index;
		result.results[r].ref_max_pos = simd_result->end_ref; // TODO: ensure these are correct with reverse
		result.results[r].query_max_pos = simd_result->end_query; // TODO: ensure these are correct with reverse
		result.results[r].score = simd_result->score;
		if (do_traceback) {
			DualAlignSIMDTraceback(simd_result, query_size, ref_size, max_tb_steps, result.results[r], start_end);

		}
		else {
			result.results[r].ref_offset = 0;
			result.results[r].query_offset = 0;
			result.results[r].total_TB_pointers = 0;
		}

		scalable_aligned_free(simd_result->trace);
		scalable_free(simd_result);

	}
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
//	std::random_shuffle(result.results.begin(), result.results.end());
}

DRAM *g_DRAM = nullptr;

InitializeScoringParameters_ptr g_InitializeScoringParameters = InitializeScoringParams;
InitializeReferenceMemory_ptr g_InitializeReferenceMemory = InitializeMemory; 
InitializeReadMemory_ptr g_InitializeReadMemory = InitializeMemory;
BatchAlignment_ptr g_BatchAlignment = BatchAlignment;
BatchAlignmentSIMD_ptr g_BatchAlignmentSIMD = BatchAlignmentSIMD;
//BatchAlignmentTwoPiece_ptr g_BatchAlignmentTwoPiece = BatchAlignmentTwoPiece;
