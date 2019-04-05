#define BOOST_LOCALE_NO_LIB
#define NOMINMAX
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>

#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <vector>
#include <algorithm>
#include <map>
#include <thread>
#include <mutex>
#include <atomic>

#include "zlib.h"
#include "kseq.h"
//#include "INIReader.h"
#include "ConfigFile.h"
#include "ntcoding.h"
#include "seed_pos_table.h"
#include "Processor.h"
#include "DRAM.h"
#include "Index.h"
#include "graph.h"

KSEQ_INIT2(, gzFile, gzread)

#define MAX_TILE_SIZE 512

std::atomic<int> num_reads(0);

//std::mutex io_lock;


float Abs(float x1) {
	float y = (x1 >= 0) ? x1 : -1 * x1;
	return y;
}

Configuration cfg;

using namespace Darwin;

struct timeval start, end_time;
long useconds, seconds, mseconds;

//LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds;
//LARGE_INTEGER Frequency;

SeedPosTable *sa;

AlignmentScoringParams params;

char* RevComp(bond::blob read) {

	size_t read_len = read.size();

	size_t read_padded = read_len;

	size_t extra = read_len % WORD_SIZE;
	if (extra != 0)
	{
		read_padded += WORD_SIZE - extra;
	}

	char* rc = (char*)scalable_aligned_malloc(read_padded, 64);

	char* seq = (char*)read.data();
	for (size_t r = 0, i = read_len; i-- > 0;) {
		if (seq[i] != 'a' && seq[i] != 'A' &&
			seq[i] != 'c' && seq[i] != 'C' &&
			seq[i] != 'g' && seq[i] != 'G' &&
			seq[i] != 't' && seq[i] != 'T' &&
			seq[i] != 'n' && seq[i] != 'N') {
			std::cerr << "Bad Nt char: " << seq[i] << std::endl;
			exit(1);
		}
		switch (seq[i]) {
		case 'a': rc[r++] = 't';
			break;

		case 'A': rc[r++] = 'T';
			break;

		case 'c': rc[r++] = 'g';
			break;

		case 'C': rc[r++] = 'G';
			break;

		case 'g': rc[r++] = 'c';
			break;

		case 'G': rc[r++] = 'C';
			break;

		case 't': rc[r++] = 'a';
			break;

		case 'T': rc[r++] = 'A';
			break;

		case 'n': rc[r++] = 'n';
			break;

		case 'N': rc[r++] = 'N';
			break;
		}
	}

	for (size_t r = read_len; r < read_padded; r++) {
		rc[r] = 'N';
	}

	return rc;
}

void PrintTileLocation(std::string& read_name, uint32_t candidate_hit, uint32_t last_hit_offset, int score, char strand) {
	size_t chr_id = std::upper_bound(Index::chr_coord.cbegin(), Index::chr_coord.cend(), candidate_hit) - Index::chr_coord.cbegin() - 1;

	std::string& chr = Index::chr_id[chr_id];

	fprintf(stderr, "%s %s %ld %ld %ld %c\n", read_name.c_str(), chr.c_str(), candidate_hit - Index::chr_coord[chr_id], last_hit_offset, score, strand);
}

void GenerateInitializeMemoryMessage(InitializeDRAMMessage& init_dram_message, const char* seq, uint32_t seq_addr, uint32_t dram_addr, int num_bytes) {
	init_dram_message.start_addr = dram_addr;
	init_dram_message.num_bytes = num_bytes;

	init_dram_message.data.resize(num_bytes / 8);

	for (int i = 0; i < num_bytes / 8; i++) {
		uint64_t data = 0;
		//memcpy(&data, seq + seq_addr + 8 * i, 8);
		for (int j = 0; j < 8; j++) {
			data = (((uint64_t)seq[seq_addr + 8 * i + j]) << 8 * j) + data;
		}
		init_dram_message.data[i] = data;
	}
}

AlignmentInputFieldsDRAM GenerateAlignmentInput(int ref_size, int query_size, int tile_size, uint64_t ref_addr, uint64_t query_addr) {

	AlignmentInputFieldsDRAM input;

	input.ref_size = ref_size;
	input.query_size = query_size;

	input.ref_bases_start_addr = ref_addr;
	input.query_bases_start_addr = query_addr;

	input.max_tb_steps = 2 * tile_size;
	input.score_threshold = 0;

	input.max_tb_steps = 2 * tile_size;
	input.score_threshold = 0;

	input.align_fields = (reverse_ref << 4) + (complement_ref << 3) + (reverse_query << 2) + (complement_query << 1) + start_end;

	return input;
}

int main(int argc, char *argv[]) {
	if (argc < 4) {
		fprintf(stderr, "Usage: testFiles.exe <REFERENCE>.fasta <READS>.fasta OVERLAP(0/1)\n");
		exit(1);
	}

    
    int do_overlap = 0;
	do_overlap = stoi(argv[3]);
//	INIReader cfg_file("params.cfg");
    ConfigFile cfg_file("params.cfg");

	fprintf(stderr, "Reading configuration ... \n");

	// GACT scoring
	cfg.gact_sub_mat[0] = cfg_file.Value("GACT_scoring", "sub_AA");
	cfg.gact_sub_mat[1] = cfg_file.Value("GACT_scoring", "sub_AC");
	cfg.gact_sub_mat[2] = cfg_file.Value("GACT_scoring", "sub_AG");
	cfg.gact_sub_mat[3] = cfg_file.Value("GACT_scoring", "sub_AT");

	cfg.gact_sub_mat[4] = cfg_file.Value("GACT_scoring", "sub_CC");
	cfg.gact_sub_mat[5] = cfg_file.Value("GACT_scoring", "sub_CG");
	cfg.gact_sub_mat[6] = cfg_file.Value("GACT_scoring", "sub_CT");

	cfg.gact_sub_mat[7] = cfg_file.Value("GACT_scoring", "sub_GG");
	cfg.gact_sub_mat[8] = cfg_file.Value("GACT_scoring", "sub_GT");

	cfg.gact_sub_mat[9] = cfg_file.Value("GACT_scoring", "sub_TT");

	cfg.gact_sub_mat[10] = cfg_file.Value("GACT_scoring", "sub_N");

	cfg.gap_open = cfg_file.Value("GACT_scoring", "gap_open");
	cfg.gap_extend = cfg_file.Value("GACT_scoring", "gap_extend");

	cfg.long_gap_open = cfg_file.Value("GACT_scoring", "long_gap_open");
	cfg.long_gap_extend = cfg_file.Value("GACT_scoring", "long_gap_extend");

	// D-SOFT parameters
	cfg.seed_size = cfg_file.Value("DSOFT_params", "seed_size");
	cfg.minimizer_window = cfg_file.Value("DSOFT_params", "minimizer_window");
	cfg.bin_size = cfg_file.Value("DSOFT_params", "bin_size");
	cfg.dsoft_threshold = cfg_file.Value("DSOFT_params", "threshold");
	cfg.num_seeds = cfg_file.Value("DSOFT_params", "num_seeds");
	cfg.seed_occurence_multiple = cfg_file.Value("DSOFT_params", "seed_occurence_multiple");
	cfg.max_candidates = cfg_file.Value("DSOFT_params", "max_candidates");
	cfg.max_stride = cfg_file.Value("DSOFT_params", "max_stride");
    cfg.do_overlap = do_overlap;
//	cfg.ignore_lower = cfg_file.GetBoolean("DSOFT_params", "ignore_lower");

	// GACT first tile
	cfg.first_tile_size = cfg_file.Value("GACT_first_tile", "first_tile_size");
	cfg.first_tile_score_threshold = cfg_file.Value("GACT_first_tile", "first_tile_score_threshold");
	cfg.first_tile_batch_size = cfg_file.Value("GACT_first_tile", "first_tile_batch_size");
	cfg.min_overlap = cfg_file.Value("GACT_first_tile", "min_overlap");
	cfg.slope_threshold = (float) cfg_file.Value("GACT_first_tile", "slope_threshold");

	// GACT extend
	cfg.tile_size = cfg_file.Value("GACT_extend", "tile_size");
	cfg.tile_overlap = cfg_file.Value("GACT_extend", "tile_overlap");
	cfg.batch_size = cfg_file.Value("GACT_extend", "batch_size");

	// Multi-threading
	cfg.num_threads = cfg_file.Value("Multithreading", "num_threads");

	// FPGA
//	cfg.processor_library = cfg_file.Get("FPGA", "processor_library");
//	cfg.num_fpgas = cfg_file.Value("FPGA", "num_fpgas");
//	cfg.chip_ids = cfg_file.Get("FPGA", "chip_ids", "0");

	//HINSTANCE hProcDLL = LoadLibraryA(cfg.processor_library.c_str());

	//size_t numProcessors;

	//if (hProcDLL != nullptr)
	//{
	//	InitializeProcessor_ptr initializeProcessor = (InitializeProcessor_ptr)GetProcAddress(hProcDLL, "InitializeProcessor");

	//	if (initializeProcessor != nullptr)
	//	{
	//		numProcessors = initializeProcessor(cfg.num_threads, cfg.num_fpgas, cfg.chip_ids);

	//			g_InitializeReferenceMemory = (InitializeReferenceMemory_ptr)GetProcAddress(hProcDLL, "InitializeReferenceMemory");
	//			g_InitializeReadMemory = (InitializeReadMemory_ptr)GetProcAddress(hProcDLL, "InitializeReadMemory");
	//			g_BatchAlignment = (BatchAlignment_ptr)GetProcAddress(hProcDLL, "BatchAlignment");
	//			//				g_BatchAlignmentTwoPiece = (BatchAlignmentTwoPiece_ptr)GetProcAddress(hProcDLL, "BatchAlignmentTwoPiece");
	//		}
	//	}
	//}

	//QueryPerformanceFrequency(&Frequency);

	//    AlignmentScoringParams params;

	// Simple scoring scheme
	params.sub_AA = cfg.gact_sub_mat[0];
	params.sub_AC = cfg.gact_sub_mat[1];
	params.sub_AG = cfg.gact_sub_mat[2];
	params.sub_AT = cfg.gact_sub_mat[3];

	params.sub_CC = cfg.gact_sub_mat[4];
	params.sub_CG = cfg.gact_sub_mat[5];
	params.sub_CT = cfg.gact_sub_mat[6];

	params.sub_GG = cfg.gact_sub_mat[7];
	params.sub_GT = cfg.gact_sub_mat[8];

	params.sub_TT = cfg.gact_sub_mat[9];

	params.sub_N = cfg.gact_sub_mat[10];

	params.gap_open = cfg.gap_open;
	params.gap_extend = cfg.gap_extend;

	params.long_gap_open = cfg.long_gap_open;
	params.long_gap_extend = cfg.long_gap_extend;

	fprintf(stderr, "Sending Alignment Parameters ...\n");
	{
		AlignmentScoringParamsResponse params_response;

		g_InitializeScoringParameters(0, params, params_response);
	}

	std::string reference_filename(argv[1]);
	std::string reads_filename(argv[2]);

	g_DRAM = new DRAM;

	Index::init();

	// Chunk size
	const size_t readBufferLimit = 1 << 6;

	{
		//// LOAD REFERENCE
		fprintf(stderr, "\nLoading reference genome ...\n");
        gettimeofday(&start, NULL);
		//QueryPerformanceCounter(&StartingTime);

		tbb::flow::graph index_graph;

		int kmer_size = cfg.seed_size;
		assert(kmer_size <= 15);
		assert(kmer_size > 3);

		size_t minimizer_window = cfg.minimizer_window;

		tbb::concurrent_vector<mini_list> minimizers;

		std::size_t histogramSize = 1ull << (kmer_size << 1);
		uint32_t* seedHistogram = (uint32_t*)scalable_calloc(histogramSize, sizeof(uint32_t));
		//#ifdef _DEBUG
		//		uint32_t* seedHistogram2 = (uint32_t*)scalable_calloc(histogramSize, sizeof(uint32_t));
		//#endif

		tbb::flow::function_node<seeder_input, seeder_input> minimizer(index_graph, tbb::flow::unlimited,
			[&](seeder_input input) {
			Read chr = get<0>(input)[0];

			uint32_t seq_len = chr.seq.size();
			char* seq_addr = (char*)chr.seq.data();
			uint32_t seq_start = seq_addr - g_DRAM->buffer;

			auto miniList = minimizers.grow_by(1);
			miniList->reserve(seq_len);

			iterate_minimizers(seq_addr, seq_len, kmer_size, minimizer_window,
				[&](uint64_t p, uint32_t m) {
				assert(m < histogramSize);

				seedHistogram[m]++;

				miniList->push_back(((uint64_t)m << 32) + p + seq_start);
			});

			//#ifdef _DEBUG
			//			mini_list miniList2;
			//			miniList2.reserve(seq_len);
			//
			//			uint32_t* window = (uint32_t*)scalable_calloc(minimizer_window, sizeof(uint32_t));
			//			uint64_t last_m = 0;
			//			uint32_t last_p = 0;
			//
			//			std::size_t rlen_2bit = (seq_len + 15) / 16;
			//			uint32_t* r_2bit = SeqToTwoBit(seq_addr, seq_len);
			//
			//			for (int p = 0; p < minimizer_window - 1; p++) {
			//				window[p] = hash32(GetSeedAtPos(r_2bit, p, kmer_size), kmer_size);
			//			}
			//
			//			for (uint32_t p = minimizer_window - 1; p < 16 * rlen_2bit - kmer_size; p++) {
			//				window[p%minimizer_window] = hash32(GetSeedAtPos(r_2bit, p, kmer_size), kmer_size);
			//				uint64_t m = Min_Window(window, minimizer_window);
			//				if ((m != last_m) || (p - last_p >= minimizer_window))
			//				{
			//					assert(m < histogramSize);
			//
			//					InterlockedIncrement(seedHistogram2 + m);
			//
			//					miniList2.push_back((m << 32) + p + seq_start);
			//					last_m = m;
			//					last_p = p;
			//				}
			//			}
			//
			//			size_t miniSize = miniList->size();
			//			assert(miniSize == miniList2.size());
			//			for (size_t i = 0; i < miniSize; i++)
			//			{
			//				assert(miniList->at(i) == miniList2[i]);
			//			}
			//
			//			scalable_free(r_2bit);
			//			scalable_free(window);
			//#endif
			return input;
		});

		tbb::flow::split_node<seeder_input, size_t> splitter(index_graph);

		tbb::flow::function_node<seeder_input, seeder_input> sender(index_graph, tbb::flow::unlimited, reference_sender_body());

		tbb::flow::make_edge(sender, splitter);

		tbb::flow::broadcast_node<seeder_input> broadcaster(index_graph);

		tbb::flow::make_edge(broadcaster, sender);

		tbb::flow::make_edge(broadcaster, minimizer);

		tbb::flow::join_node<seeder_input> gatekeeper(index_graph);

		tbb::flow::make_edge(gatekeeper, broadcaster);

		tbb::flow::buffer_node<size_t> ticketer(index_graph);

		// Allocate tickets
		for (size_t t = 0ull; t < cfg.num_threads; t++)
			//for (size_t t = 0ull; t < 1ull; t++)
			ticketer.try_put(t);

		tbb::flow::make_edge(tbb::flow::output_port<1>(splitter), ticketer);

		tbb::flow::make_edge(ticketer, tbb::flow::input_port<1>(gatekeeper));

		gzFile f_rd = gzopen(argv[1], "r");
		if (!f_rd) { fprintf(stderr, "cant open file: %s\n", argv[1]); exit(EXIT_FAILURE); }

		kseq_t *kseq_rd = kseq_init(f_rd);

		tbb::flow::source_node<reader_output> reader(index_graph,
			[&](reader_output &reads) -> bool {

			reads.clear();

			if (kseq_read(kseq_rd) >= 0)
			{
				size_t seq_len = kseq_rd->seq.l;
                size_t seq_len_unpadded = seq_len;

				if (seq_len > readBufferLimit) {

					memcpy(g_DRAM->buffer + g_DRAM->referenceSize, kseq_rd->seq.s, seq_len);

					Read read;

					read.description = std::string(kseq_rd->name.s, kseq_rd->name.l);

					read.seq.assign(g_DRAM->buffer + g_DRAM->referenceSize, seq_len);

					size_t extra = seq_len % WORD_SIZE;
					if (extra != 0)
					{
						extra = WORD_SIZE - extra;
						memset(g_DRAM->buffer + g_DRAM->referenceSize + seq_len, 'N', extra);

						seq_len += extra;
					}

					reads.push_back(read);

					g_DRAM->referenceSize += seq_len;

					Index::chr_id.push_back(read.description);

					Index::chr_len.push_back(seq_len);
                    Index::chr_len_unpadded.push_back(seq_len_unpadded);

					Index::chr_coord.push_back(g_DRAM->referenceSize);

					return true;
				}
				else {
					return false;
				}
			}

			return false;
		});

		tbb::flow::make_edge(reader, tbb::flow::input_port<0>(gatekeeper));

		index_graph.wait_for_all();

		//#ifdef _DEBUG
		//		for (size_t i = 0; i < histogramSize; i++)
		//		{
		//			assert(seedHistogram[i] == seedHistogram2[i]);
		//		}
		//
		//		scalable_free(seedHistogram2);
		//#endif

		gzclose(f_rd);

		g_DRAM->bufferPosition = g_DRAM->referenceSize;

		fprintf(stderr, "Reference length: %lld\n", g_DRAM->referenceSize);
        
        gettimeofday(&end_time, NULL);

        useconds = end_time.tv_usec - start.tv_usec;
        seconds = end_time.tv_sec - start.tv_sec;
        mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;

        std::cerr << "Time elapsed (loading reference): " << mseconds <<" msec" << std::endl;


		//QueryPerformanceCounter(&EndingTime);

		//ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
		//ElapsedMicroseconds.QuadPart *= 1000000;
		//ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
//		fprintf(stderr, "Time elapsed (loading reference genome): %lld msec\n", (ElapsedMicroseconds.QuadPart / 1000));

		//// FINALIZE SEED POSITION TABLE
		fprintf(stderr, "\nFinalizing seed position table ...\n");
        gettimeofday(&start, NULL);
		//QueryPerformanceCounter(&StartingTime);

		sa = new SeedPosTable(g_DRAM->referenceSize, cfg.seed_size, cfg.minimizer_window, cfg.max_stride, cfg.seed_occurence_multiple, cfg.bin_size,
			minimizers, seedHistogram, histogramSize);

		scalable_free(seedHistogram);
        
        gettimeofday(&end_time, NULL);

        useconds = end_time.tv_usec - start.tv_usec;
        seconds = end_time.tv_sec - start.tv_sec;
        mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;

        std::cerr << "Time elapsed (finalizing seed position table): " << mseconds <<" msec" << std::endl;


		//QueryPerformanceCounter(&EndingTime);

		//ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
		//ElapsedMicroseconds.QuadPart *= 1000000;
		//ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
		//fprintf(stderr, "Time elapsed (seed position table construction): %lld msec\n", (ElapsedMicroseconds.QuadPart / 1000));
	}

	//// SEND REFERENCE
	//QueryPerformanceCounter(&StartingTime);
//	fprintf(stderr, "\nSending reference ...\n");

	//int word_size = WORD_SIZE;
	//int num_bytes_to_send = 0;
	//size_t max_char_to_send = MAX_CHAR_TO_SEND;

	//int num_words = 0;
	//uint32_t dram_start_addr = 0;
	//uint32_t reference_start_addr = 0;
	//uint32_t reads_start_addr = 0;

	////TODO
	//fpga_writer_lock.lock();

	//for (dram_start_addr = 0; dram_start_addr < g_DRAM->referenceSize; dram_start_addr += num_bytes_to_send) {
	//	// since reference sequences are already padded to cfg.bin_size 
	//	// reference should be padded to word_size
	//	num_bytes_to_send = std::min(max_char_to_send, g_DRAM->referenceSize - dram_start_addr);

	//	assert(num_bytes_to_send >= word_size);

	//	InitializeDRAMMessage   init_dram_message;
	//	GenerateInitializeMemoryMessage(init_dram_message, g_DRAM->buffer, dram_start_addr, dram_start_addr, num_bytes_to_send);
	//	//            WriteHaasMessage ("InitializeDRAMRequest." + to_string(num_words++) + ".bin", Function::InitializeDRAM, Magic::InitializeDRAM, init_dram_message);

	//	{
	//		InitializeDRAMMessageResponse init_dram_response;
	//		g_InitializeMemory(0, g_DRAM->buffer, init_dram_message, init_dram_response);
	//	}
	//	//            dram_start_addr += word_size;
	//}

	//fpga_writer_lock.unlock();

	//QueryPerformanceCounter(&EndingTime);

	//ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
	//ElapsedMicroseconds.QuadPart *= 1000000;
	//ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
	//fprintf(stderr, "Time elapsed (sending reference): %lld msec\n", (ElapsedMicroseconds.QuadPart / 1000));

	//// CONSTRUCT SEED POSITION TABLE
//	fprintf(stderr, "\nConstructing seed position table ...\n");
	//QueryPerformanceCounter(&StartingTime);

	//sa = new SeedPosTable(g_DRAM->buffer, g_DRAM->referenceSize, cfg.seed_size, cfg.minimizer_window, cfg.seed_occurence_multiple, cfg.bin_size);

	//QueryPerformanceCounter(&EndingTime);

	//ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
	//ElapsedMicroseconds.QuadPart *= 1000000;
	//ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
//	fprintf(stderr, "Time elapsed (seed position table construction): %lld msec\n", (ElapsedMicroseconds.QuadPart / 1000));

	fprintf(stderr, "\nAligning reads ...\n");
    gettimeofday(&start, NULL);
	//QueryPerformanceCounter(&StartingTime);

	tbb::flow::graph align_graph;

	tbb::flow::function_node<printer_input, size_t> printer(align_graph, tbb::flow::unlimited, printer_body());

	extender_node extender(align_graph, tbb::flow::unlimited, extender_body());

	tbb::flow::make_edge(tbb::flow::output_port<0>(extender), printer);
//	tbb::flow::make_edge(extender, printer);

	tbb::flow::function_node<filter_input, extender_input> filter(align_graph, tbb::flow::unlimited, filter_body());

	tbb::flow::make_edge(filter, extender);

	tbb::flow::function_node<seeder_input, filter_input> seeder(align_graph, tbb::flow::unlimited, seeder_body());

	tbb::flow::make_edge(seeder, filter);

	tbb::flow::function_node<seeder_input, seeder_input> sender(align_graph, tbb::flow::unlimited, reference_sender_body());

	tbb::flow::make_edge(sender, seeder);

	tbb::flow::join_node<seeder_input> gatekeeper(align_graph);

	tbb::flow::make_edge(gatekeeper, sender);

	tbb::flow::buffer_node<size_t> ticketer(align_graph);

	// Allocate tickets
	for (size_t t = 0ull; t < cfg.num_threads; t++)
		//for (size_t t = 0ull; t < 1ull; t++)
		ticketer.try_put(t);

	tbb::flow::make_edge(tbb::flow::output_port<1>(extender), ticketer);

	tbb::flow::make_edge(ticketer, tbb::flow::input_port<1>(gatekeeper));

	gzFile f_rd = gzopen(argv[2], "r");
	if (!f_rd) { fprintf(stderr, "cant open reads file: %s\n", argv[2]); exit(EXIT_FAILURE); }

	kseq_t *kseq_rd = kseq_init(f_rd);

	tbb::flow::source_node<reader_output> reader(align_graph,
		[&](reader_output &reads) -> bool {

		size_t readBufferSize = 0;

		reads.clear();

		while (true)
		{
			if (kseq_read(kseq_rd) >= 0)
			{
				num_reads += 1;

				// align on word size
				size_t extra = g_DRAM->bufferPosition % WORD_SIZE;
				if (extra != 0)
				{
					extra = WORD_SIZE - extra;

					g_DRAM->bufferPosition += extra;
				}

				size_t seq_len = kseq_rd->seq.l;

				if (seq_len > readBufferLimit) {
					// Wrap around if we have reached the end of the buffer
					if (g_DRAM->bufferPosition + WORD_SIZE + seq_len > g_DRAM->size)
					{
						g_DRAM->bufferPosition = g_DRAM->referenceSize;
					}

					std::memcpy(g_DRAM->buffer + g_DRAM->bufferPosition, kseq_rd->seq.s, seq_len);

					Read read;

					read.description = std::string(kseq_rd->name.s, kseq_rd->name.l);

					read.seq = bond::blob(g_DRAM->buffer + g_DRAM->bufferPosition, seq_len);
					char *rev_read_char = RevComp(read.seq);
					read.rc_seq = bond::blob(rev_read_char, seq_len);

					extra = seq_len % WORD_SIZE;
					if (extra != 0)
					{
						extra = WORD_SIZE - extra;
						memset(g_DRAM->buffer + g_DRAM->bufferPosition + seq_len, 'N', extra);

						seq_len += extra;
					}

					g_DRAM->bufferPosition += seq_len;

					readBufferSize += seq_len;

					reads.push_back(read);
				}

				if ((readBufferSize > readBufferLimit) || (kseq_rd->f->is_eof))
				{
					return true;
				}
			}
			else
			{
				return false;
			}
		}
	}, true);

	tbb::flow::make_edge(reader, tbb::flow::input_port<0>(gatekeeper));

	align_graph.wait_for_all();

	gzclose(f_rd);

	//QueryPerformanceCounter(&EndingTime);

	//ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
	//ElapsedMicroseconds.QuadPart *= 1000000;
	//ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
	//fprintf(stderr, "Time elapsed (aligning reads): %lld msec\n\n", (ElapsedMicroseconds.QuadPart / 1000));

	std::cerr << "#reads: " << num_reads << std::endl;
	std::cerr << "#filter tiles: " << filter_body::num_filter_tiles << std::endl;
	std::cerr << "#extend requests: " << filter_body::num_extend_requests << std::endl;
	std::cerr << "#slope filtered: " << filter_body::num_slope_filtered << std::endl;
	std::cerr << "#extend tiles: " << extender_body::num_extend_tiles << std::endl;
	std::cerr << "#active tiles: " << extender_body::num_active_tiles << std::endl;
	std::cerr << "#large tiles: " << extender_body::num_large_tiles << std::endl;

    gettimeofday(&end_time, NULL);

    useconds = end_time.tv_usec - start.tv_usec;
    seconds = end_time.tv_sec - start.tv_sec;
    mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    
    std::cerr << "Time elapsed (aligning reads): " << mseconds <<" msec" << std::endl;

	return 0;
}

