#pragma once
#include "Darwin_reflection.h"

#include "DRAM.h"
#include "Index.h"
#include <tbb/scalable_allocator.h>
#include "immintrin.h"

#define INF (1 << 30)
#define MAX_TILE_SIZE 512

typedef int AlnOp;
enum AlnOperands { ZERO_OP, INSERT_OP, DELETE_OP, LONG_INSERT_OP, LONG_DELETE_OP};
enum states { Z, I, D, M , L_Z, L_I, L_D};

/* for traceback */
#define TRACEBACK_T_MASK     8160  /* all bits set except the first five */
#define TRACEBACK_F_MASK     7807  /* all bits set except the F bits */
#define TRACEBACK_F_L_MASK   2047  /* all bits set except the F_L bits */

#define TRACEBACK_ZERO_T        0
#define TRACEBACK_DEL_T         1
#define TRACEBACK_INS_T         2
#define TRACEBACK_DEL_L_T       4
#define TRACEBACK_INS_L_T       8
#define TRACEBACK_DIAG          16
#define TRACEBACK_DIAG_DEL      32
#define TRACEBACK_DEL           64
#define TRACEBACK_DIAG_INS      128
#define TRACEBACK_INS           256
#define TRACEBACK_DIAG_DEL_L    512
#define TRACEBACK_DEL_L         1024
#define TRACEBACK_DIAG_INS_L    2048
#define TRACEBACK_INS_L         4096

typedef union __m256i_16 {
    __m256i m;
    int16_t v[16];
} __m256i_16_t;

typedef struct result {
    int score;      /* alignment score */
    int end_query;  /* end position of query sequence */
    int end_ref;    /* end position of reference sequence */
    __m256i *trace;
} result_t;

using namespace Darwin;

typedef size_t(*InitializeProcessor_ptr)(int t, int f, std::string s);
typedef void(*InitializeScoringParameters_ptr)(size_t token, AlignmentScoringParams& request, AlignmentScoringParamsResponse& response);
typedef void(*InitializeReferenceMemory_ptr)(size_t token, char* dram, InitializeDRAMMessage& request, InitializeDRAMMessageResponse& response);
typedef void(*InitializeReadMemory_ptr)(size_t token, char* dram, InitializeDRAMMessage& request, InitializeDRAMMessageResponse& response);
typedef void(*BatchAlignment_ptr)(size_t token, char* dram, BatchAlignmentInputFieldsDRAM& request, BatchAlignmentResultDRAM& result);
typedef void(*BatchAlignmentSIMD_ptr)(size_t token, char* dram, BatchAlignmentInputFieldsDRAM& request, BatchAlignmentResultDRAM& result);

extern DRAM *g_DRAM;

extern InitializeScoringParameters_ptr g_InitializeScoringParameters;
extern InitializeReferenceMemory_ptr g_InitializeReferenceMemory;
extern InitializeReadMemory_ptr g_InitializeReadMemory;
extern BatchAlignment_ptr g_BatchAlignment;
extern BatchAlignmentSIMD_ptr g_BatchAlignmentSIMD;
//extern BatchAlignmentTwoPiece_ptr g_BatchAlignmentTwoPiece;
