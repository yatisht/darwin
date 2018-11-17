#pragma once
#include "Darwin_reflection.h"

#include "DRAM.h"
#include "Index.h"

#define INF (1 << 24)
#define MAX_TILE_SIZE 512

typedef int AlnOp;
enum AlnOperands { ZERO_OP, INSERT_OP, DELETE_OP, LONG_INSERT_OP, LONG_DELETE_OP};
enum states { Z, I, D, M , L_Z, L_I, L_D};

using namespace Darwin;

typedef size_t(*InitializeProcessor_ptr)(int t, int f, std::string s);
typedef void(*InitializeScoringParameters_ptr)(size_t token, AlignmentScoringParams& request, AlignmentScoringParamsResponse& response);
typedef void(*InitializeReferenceMemory_ptr)(size_t token, char* dram, InitializeDRAMMessage& request, InitializeDRAMMessageResponse& response);
typedef void(*InitializeReadMemory_ptr)(size_t token, char* dram, InitializeDRAMMessage& request, InitializeDRAMMessageResponse& response);
typedef void(*BatchAlignment_ptr)(size_t token, char* dram, BatchAlignmentInputFieldsDRAM& request, BatchAlignmentResultDRAM& result);
//typedef void(*BatchAlignmentTwoPiece_ptr)(size_t token, char* dram, BatchAlignmentInputFieldsDRAM& request, BatchAlignmentResultDRAM& result);
//typedef void(*BatchAlignmentTwoPiece_ptr)(size_t token, char* dram, BatchAlignmentInputFieldsDRAM& request, BatchAlignmentTwoPieceResultDRAM& result);

extern DRAM *g_DRAM;

extern InitializeScoringParameters_ptr g_InitializeScoringParameters;
extern InitializeReferenceMemory_ptr g_InitializeReferenceMemory;
extern InitializeReadMemory_ptr g_InitializeReadMemory;
extern BatchAlignment_ptr g_BatchAlignment;
//extern BatchAlignmentTwoPiece_ptr g_BatchAlignmentTwoPiece;
