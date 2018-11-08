#include "graph.h"
std::mutex io_lock;

int mat_offset[] = { 0, 1, 3, 6 };

extern int do_overlap;

Read maf_printer_body::operator()(printer_input input)
{
    Read read = input.first;
    ExtendAlignments e = input.second;

    int score = 0;
    int index = 0;
    int open = 0;
    int short_gap_penalty = 0;
    int long_gap_penalty = 0;
    for (int l = 0; l < e.aligned_reference_str.length(); l++) {
        char r = e.aligned_reference_str[l];
        char q = e.aligned_query_str[l];
        if ((r == '-') || (q == '-')) {
            short_gap_penalty += (open) ? cfg.gap_extend : cfg.gap_open;
            long_gap_penalty += (open) ? cfg.long_gap_extend : cfg.long_gap_open;
            open = 1;
        }
        else {
            int r_nt = NtChar2Int(r);
            int q_nt = NtChar2Int(q);
            if ((r_nt <= 3) && (q_nt <= 3)) {
                if (r_nt > q_nt) {
                    index = q_nt * 4 + r_nt - mat_offset[q_nt];
                }
                else {
                    index = r_nt * 4 + q_nt - mat_offset[r_nt];
                }
                score += cfg.gact_sub_mat[index];
            }
            else {
                score += cfg.gact_sub_mat[10];
            }
            score += (long_gap_penalty < short_gap_penalty) ? short_gap_penalty : long_gap_penalty;
            open = 0;
            short_gap_penalty = 0;
            long_gap_penalty = 0;
        }
    }
    uint32_t ref_start = 1 + e.reference_start_offset;
    uint32_t query_start = 1 + e.query_start_offset;
    uint32_t ref_align_len = ((e.reference_end_offset + 1) - e.reference_start_offset);
    uint32_t query_align_len = ((e.query_end_offset + 1) - e.query_start_offset);

		io_lock.lock();
		// Alignment spans 80% of the read
		bool do_print = (do_overlap) ? ((5 * query_align_len > 4 * read.seq.size()) || (ref_align_len >= cfg.min_overlap) || (query_align_len >= cfg.min_overlap)) : (5 * query_align_len > 4 * read.seq.size());
		if (do_print) {
			printf("a score=%d\n", score);
			printf("s\t%s\t%lu\t%lu\t+\t%lu\t", Index::chr_id[e.chr_id].c_str(), ref_start, ref_align_len, Index::chr_len[e.chr_id]);
			puts(e.aligned_reference_str.c_str());
			printf("\n");
			printf("s\t%s\t%lu\t%lu\t%c\t%lu\t", read.description.c_str(), query_start, query_align_len, e.strand, read.seq.size());
			puts(e.aligned_query_str.c_str());
			printf("\n\n");
		}
		io_lock.unlock();

    //assert(score >= cfg.first_tile_score_threshold);

    return read;
};
