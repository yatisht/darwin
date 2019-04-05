#pragma once
#define BOOST_LOCALE_NO_LIB
#define NOMINMAX
#include <algorithm>

#include <tbb/flow_graph.h>
#include <tbb/reader_writer_lock.h>
#include <tbb/scalable_allocator.h>
#include <mutex>

#include "ntcoding.h"
#include "seed_pos_table.h"
#include "Processor.h"
#include "DRAM.h"
#include "Index.h"

#define TB_MASK 3

#define MAX_CHAR_TO_SEND 2048

//#define traceback 0
#define reverse_ref (1 << 4)
#define complement_ref (1 << 3)
#define reverse_query (1 << 2)
#define complement_query (1 << 1)
#define start_end 1

struct Configuration {
	// GACT scoring
	int gact_sub_mat[11];
	int gap_open;
	int gap_extend;
	int long_gap_open;
	int long_gap_extend;

	// D-SOFT parameters
	int seed_size;
	int minimizer_window;
	uint32_t bin_size;
	int dsoft_threshold;
	int num_seeds;
	int seed_occurence_multiple;
	int max_candidates;
	int num_nz_bins;
    int max_stride;
    int do_overlap;
	bool ignore_lower;

	// GACT first tile
	int first_tile_size;
	int first_tile_score_threshold;
	int first_tile_batch_size;

	int min_overlap;
	float slope_threshold;

	//GACT extend 
	int tile_size;
	int tile_overlap;
	int batch_size;

	//Multi-threading
	int num_threads;

    //FPGA
	std::string processor_library;
    int num_fpgas;
//    std::string chip_ids;
};

extern Configuration cfg;

extern SeedPosTable *sa;

void GenerateInitializeMemoryMessage(InitializeDRAMMessage& init_dram_message, const char* seq, uint32_t seq_addr, uint32_t dram_addr, int num_bytes);

struct Read {
	std::string description;
	bond::blob seq;
	bond::blob rc_seq;
};

struct ExtendLocations {
	int read_num;
	int chr_id;
	int score;
	uint32_t reference_pos;
	uint32_t query_pos;
    vector<uint64_t> left_hit_offsets;
	vector<uint64_t> right_hit_offsets;
};

static inline bool CompareExtendLocations(ExtendLocations e1, ExtendLocations e2) {
	return e1.score > e2.score;
}

struct ExtendAlignments {
	int read_num;
	int chr_id;
	uint32_t curr_reference_offset;
	uint32_t curr_query_offset;
	uint32_t reference_start_offset;
	uint32_t query_start_offset;
	uint32_t reference_end_offset;
	uint32_t query_end_offset;
	uint32_t reference_start_addr;
	uint32_t query_start_addr;
	uint32_t reference_length;
	uint32_t query_length;
	int left_extension_done;
	int right_extension_done;
    bool used_large_tile;
    bool do_print;
	char strand;
	std::string aligned_reference_str;
	std::string aligned_query_str;
    int score;
    int chain_score;
    vector<uint64_t> left_hit_offsets;
	vector<uint64_t> right_hit_offsets;
};

struct SamFields {
    std::string QNAME;
    uint16_t FLAG;
    std::string RNAME;
    size_t  POS;
    uint8_t MAPQ;
    std::string CIGAR;
    std::string RNEXT;
    size_t PNEXT;
    int32_t TLEN;
    std::string SEQ;
    std::string QUAL;
    int32_t SCORE;
    int32_t CHAIN_SCORE;
};

typedef std::vector<Read> reader_output;

typedef tbb::flow::tuple<reader_output, size_t> seeder_input;

struct seeder_data
{
	std::vector<std::size_t> fwAnchorBuckets;
	std::vector<std::size_t> rcAnchorBuckets;

	std::vector<Anchors> fwAnchors;
	std::vector<Anchors> rcAnchors;
};

typedef tbb::flow::tuple<seeder_data, size_t> seeder_output;

typedef tbb::flow::tuple<reader_output, seeder_data> filter_payload;
typedef tbb::flow::tuple<filter_payload, size_t> filter_input;

struct filter_data
{
	std::vector<ExtendLocations> fwLocations;
	std::vector<ExtendLocations> rcLocations;
};

typedef tbb::flow::tuple<reader_output, filter_data> extender_payload;
typedef tbb::flow::tuple<extender_payload, size_t> extender_input;

struct extend_data
{
	std::vector<ExtendAlignments> extend_alignments;
};
typedef tbb::flow::tuple<reader_output, extend_data> printer_payload;
typedef tbb::flow::tuple<printer_payload, size_t> printer_input;
//typedef std::pair<Read, ExtendAlignments> printer_input;
//typedef std::pair<Read, size_t> printer_input;

typedef tbb::flow::tuple<printer_input, size_t> extender_output;

typedef tbb::flow::multifunction_node<extender_input, extender_output> extender_node;

struct printer_body
{
	static std::atomic<int> done_header;
    SamFields AlignmentToSam (Read r, ExtendAlignments e);
    void sam_printer(printer_input input);
    void mhap_printer(printer_input input);
	size_t operator()(printer_input input);
};

//extern tbb::reader_writer_lock fpga_writer_lock;

extern std::mutex io_lock;

struct reference_sender_body
{
	seeder_input operator()(seeder_input input);
};

struct read_sender_body
{
	seeder_input operator()(seeder_input input);
};


struct seeder_body
{
	filter_input operator()(seeder_input input);
};

struct filter_body
{
	static std::atomic<int> num_filter_tiles;
	static std::atomic<int> num_extend_requests;
	static std::atomic<int> num_slope_filtered;

	extender_input operator()(filter_input input);
	void slopeFilter(std::deque<ExtendLocations> &read_extend_locations, std::vector<ExtendLocations> &extend_locations);
};

struct extender_body
{
	static std::atomic<int> num_extend_tiles;
	static std::atomic<int> num_active_tiles;
	static std::atomic<int> num_large_tiles;

//	printer_input operator()(extender_input input, extender_node::output_ports_type &op);
	void operator()(extender_input input, extender_node::output_ports_type &op);
	ExtendAlignments makeForwardAlignment(std::vector<Read> &batch, std::vector<ExtendLocations>::const_iterator &loc);
	ExtendAlignments makeBackwardAlignment(std::vector<Read> &batch, std::vector<ExtendLocations>::const_iterator &loc);
    int AlignmentScore(std::string r, std::string q);
};
