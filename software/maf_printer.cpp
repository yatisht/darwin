#include "graph.h"
std::mutex io_lock;

extern int do_overlap;

size_t maf_printer_body::operator()(printer_input input)
{
    auto &payload = get<0>(input); 

    auto &reads = get<0>(payload);

    auto &data = get<1>(payload);

    size_t token = get<1>(input);
    //
    Read read;

    auto& extend_alignments = data.extend_alignments;
    std::stable_sort(extend_alignments.begin(), extend_alignments.end(),
            [](const ExtendAlignments &e1, const ExtendAlignments &e2) {
            return ((e1.read_num < e2.read_num) || ((e1.read_num == e2.read_num) && (e1.score > e2.score)));
            });

    for (auto e1 = extend_alignments.begin(); e1 != extend_alignments.end(); e1++) {
        if (!e1->do_print) {
            continue;
        }
        uint32_t s_1 = e1->query_start_offset;
        uint32_t e_1 = e1->query_end_offset;
        for (auto e2 = e1 + 1; e2 != extend_alignments.end(); e2++) {
            if (!e2->do_print) {
                continue;
            }
            if (e2->read_num != e1->read_num) {
                break;
            }
            uint32_t s_2 = e2->query_start_offset;
            uint32_t e_2 = e2->query_end_offset;
            uint32_t s = std::max(s_1, s_2);
            uint32_t e = std::min(e_1, e_2);
            uint32_t overlap = 0;
            if (e > s) {
                overlap = e-s;
            }
            if (2*overlap > (e_2-s_2)) {
                e2->do_print = false;
            }
        }
    }

    for (auto e: extend_alignments) {
        if (e.do_print) {
            Read read = reads[e.read_num];
            int score = e.score;
            uint32_t ref_start = 1 + e.reference_start_offset;
            uint32_t query_start = 1 + e.query_start_offset;
            uint32_t ref_align_len = ((e.reference_end_offset + 1) - e.reference_start_offset);
            uint32_t query_align_len = ((e.query_end_offset + 1) - e.query_start_offset);

            io_lock.lock();
            // Alignment spans 80% of the read
//            bool do_print = (do_overlap) ? ((5 * query_align_len > 4 * read.seq.size()) || (ref_align_len >= cfg.min_overlap) || (query_align_len >= cfg.min_overlap)) : (5 * query_align_len > 4 * read.seq.size());
//            if (do_print) {
                printf("a score=%d\n", score);
                printf("s\t%s\t%lu\t%lu\t+\t%lu\t", Index::chr_id[e.chr_id].c_str(), ref_start, ref_align_len, Index::chr_len[e.chr_id]);
                puts(e.aligned_reference_str.c_str());
                printf("s\t%s\t%lu\t%lu\t%c\t%lu\t", read.description.c_str(), query_start, query_align_len, e.strand, read.seq.size());
                puts(e.aligned_query_str.c_str());
                printf("\n");
//            }
            io_lock.unlock();

        }
    }

    return token;
};
