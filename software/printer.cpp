#include "graph.h"
#include <atomic>
std::mutex io_lock;

std::atomic<int> printer_body::done_header(0);

void printer_body::sam_printer(printer_input input) {
    auto &payload = get<0>(input); 

    auto &reads = get<0>(payload);

    auto &data = get<1>(payload);

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
            // SAM Format
            SamFields sf = AlignmentToSam(read, e);
            
            io_lock.lock();
            if (done_header == 0) {
                std::cout << "@HD\tVN:1.6\tSO:coordinate\n";
                int num_chrom = Index::chr_id.size();
                for (int c = 0; c < num_chrom; c++) {
                    std::cout << "@SQ\tSN:" << Index::chr_id[c] << "\tLN:" << Index::chr_len_unpadded[c] << "\n"; 
                }
                done_header = 1;
            }


            std::cout << sf.QNAME << "\t";
            std::cout << sf.FLAG << "\t";
            std::cout << sf.RNAME << "\t";
            std::cout << sf.POS << "\t";
            std::cout << (int)sf.MAPQ << "\t";
            std::cout << sf.CIGAR << "\t";
            std::cout << sf.RNEXT << "\t";
            std::cout << sf.PNEXT << "\t";
            std::cout << sf.TLEN << "\t";
            std::cout << sf.SEQ << "\t";
            std::cout << sf.QUAL << "\t";
            std::cout << "AS:i:" << sf.SCORE << "\t";
            std::cout << "ZS:i:" << sf.CHAIN_SCORE;
            std::cout << "\n";

            //MAF format
            //            int score = e.score;
            //            uint32_t ref_start = 1 + e.reference_start_offset;
            //            uint32_t query_start = 1 + e.query_start_offset;
            //            uint32_t ref_align_len = ((e.reference_end_offset + 1) - e.reference_start_offset);
            //            uint32_t query_align_len = ((e.query_end_offset + 1) - e.query_start_offset);
            //            printf("a score=%d\n", score);
            //            printf("s\t%s\t%lu\t%lu\t+\t%lu\t", Index::chr_id[e.chr_id].c_str(), ref_start, ref_align_len, Index::chr_len[e.chr_id]);
            //            puts(e.aligned_reference_str.c_str());
            //            printf("s\t%s\t%lu\t%lu\t%c\t%lu\t", read.description.c_str(), query_start, query_align_len, e.strand, read.seq.size());
            //            puts(e.aligned_query_str.c_str());
            //            printf("\n");

            io_lock.unlock();
        }
    }
}

void printer_body::mhap_printer(printer_input input) {
    auto &payload = get<0>(input); 

    auto &reads = get<0>(payload);

    auto &data = get<1>(payload);

    //
    Read read;

    auto& extend_alignments = data.extend_alignments;
    std::stable_sort(extend_alignments.begin(), extend_alignments.end(),
            [](const ExtendAlignments &e1, const ExtendAlignments &e2) {
            return ((e1.read_num < e2.read_num) || ((e1.read_num == e2.read_num) && (e1.chr_id < e2.chr_id)) || ((e1.read_num == e2.read_num) && (e1.chr_id == e2.chr_id) && (e1.score > e2.score)));
            });

    for (auto e1 = extend_alignments.begin(); e1 != extend_alignments.end(); e1++) {
        
        uint32_t ref_end = 1 + e1->reference_end_offset;
        uint32_t query_end = 1 + e1->query_end_offset;
        uint32_t ref_len = e1->reference_length; 
        uint32_t query_len = e1->query_length; 
        if ((ref_end < (9*ref_len)/10) && (query_end < (9*query_len)/10)) {
            e1->do_print = false;
        }
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
            if (e1->chr_id != e2->chr_id) {
                break;
            }
            e2->do_print = false;
        }
    }

    for (auto e: extend_alignments) {
        if (e.do_print) {
            Read read = reads[e.read_num];
            io_lock.lock();
            
            int score = e.score;
            std::string r1_name = Index::chr_id[e.chr_id];
            std::string r2_name = read.description; 
            int r2_strand = (e.strand == '-') ? 1 : 0;
            uint32_t ref_start = 1 + e.reference_start_offset;
            uint32_t query_start = 1 + e.query_start_offset;
            uint32_t ref_end = 1 + e.reference_end_offset;
            uint32_t query_end = 1 + e.query_end_offset;
            uint32_t ref_align_len = ((e.reference_end_offset + 1) - e.reference_start_offset);
            uint32_t query_align_len = ((e.query_end_offset + 1) - e.query_start_offset);
            uint32_t ovl_len = 0, matches = 0;
            for (int i  = 0; i < e.aligned_reference_str.length(); i++) {
                if (toupper(e.aligned_reference_str[i]) == toupper(e.aligned_query_str[i])) {
                    matches +=1;
                }
            }
            ovl_len = (ref_align_len + query_align_len)/2;
            float error = (1.0*(ovl_len - matches))/ovl_len;

            if ((ovl_len >= cfg.min_overlap) && (r1_name != r2_name)){
                printf ("%s %s %.3f %d %d %lu %lu %lu %d %lu %lu %lu\n", r1_name.c_str(), r2_name.c_str(), error, matches, 0, ref_start, ref_end, Index::chr_len_unpadded[e.chr_id], r2_strand, query_start, query_end, read.seq.size());
                printf("%s\n", e.aligned_reference_str.c_str());
                printf("%s\n", e.aligned_query_str.c_str());
                printf ("%s %s %.3f %d %d %lu %lu %lu %d %lu %lu %lu\n", r2_name.c_str(), r1_name.c_str(), error, matches, r2_strand, query_start, query_end, read.seq.size(), 0, ref_start, ref_end, Index::chr_len_unpadded[e.chr_id]);
                printf("%s\n", e.aligned_query_str.c_str());
                printf("%s\n", e.aligned_reference_str.c_str());

            }
            io_lock.unlock();
        }
    }
}
size_t printer_body::operator()(printer_input input)
{
    size_t token = get<1>(input);
    if (cfg.do_overlap == 1) {
        mhap_printer(input);
    }
    else {
        sam_printer(input);
    }
    return token;
};


SamFields printer_body::AlignmentToSam (Read read, ExtendAlignments e) {
    SamFields sf;

    sf.QNAME  = std::string(read.description.c_str());
    sf.FLAG   = 0;
    if (e.strand == '-') {
        sf.FLAG += 16;
        char* read_char = (char*) read.rc_seq.data(); 
        const size_t read_len = read.rc_seq.size(); 
        sf.SEQ = std::string(read_char, read_len);
    }
    else {
        char* read_char = (char*) read.seq.data();
        const size_t read_len = read.seq.size(); 
        sf.SEQ = std::string(read_char, read_len);
    }

    sf.FLAG += 64;

    sf.RNAME = std::string(Index::chr_id[e.chr_id].c_str()); 
    sf.POS = 1 + e.reference_start_offset;

    // UNUNSED. FOR LATER
    sf.MAPQ = 60;

    // ------------ CIGAR ------------- //
    std::vector<char> cigar;
    cigar.clear();

    char op = 'Z';
    char prev_op = 'Z'; 
    size_t num_op = 0;
    std::string num_op_str = "";

    if (e.query_start_offset > 0) {
        num_op = e.query_start_offset;
        op = 'S';
        num_op_str = std::to_string(num_op);
        for (int i = 0; i < num_op_str.length(); i++) {
            cigar.push_back(num_op_str[i]);
        }
        cigar.push_back(op);
    }

    assert(e.aligned_reference_str.length() == e.aligned_query_str.length());
    num_op = 0;
    for (int i  = 0; i < e.aligned_reference_str.length(); i++) {
        char r = e.aligned_reference_str[i];
        char q = e.aligned_query_str[i];
        if (r == '-') {
            op = 'I';
        }
        else if (q == '-') {
            op = 'D';
        }
        else {
            op = 'M';
        }
//        else if (toupper(r) == toupper(q)) {
//            op = '=';
//        }
//        else {
//            op = 'X';
//        }
        if (op == prev_op) {
            num_op++;
        }
        else {
            if (num_op > 0) {
                num_op_str = std::to_string(num_op);
                for (int i = 0; i < num_op_str.length(); i++) {
                    cigar.push_back(num_op_str[i]);
                }
                cigar.push_back(prev_op);
            }
            num_op = 1;
        }
        prev_op = op;
    }

    if (num_op > 0) {
        num_op_str = std::to_string(num_op);
        for (int i = 0; i < num_op_str.length(); i++) {
            cigar.push_back(num_op_str[i]);
        }
        cigar.push_back(prev_op);
    }

    num_op = (e.query_length - e.query_end_offset - 1);
    if (num_op > 0) {
        op = 'S';
        num_op_str = std::to_string(num_op);
        for (int i = 0; i < num_op_str.length(); i++) {
            cigar.push_back(num_op_str[i]);
        }
        cigar.push_back(op);
    }

    sf.CIGAR = (cigar.size() == 0) ? "*" : std::string(cigar.begin(), cigar.end());


    // ------------ END CIGAR ------------- //

    sf.RNEXT = "*";
    sf.PNEXT = 0;
    sf.TLEN = 0;
    sf.QUAL = "*";
    sf.SCORE = e.score;
    //FOR LATER
    sf.CHAIN_SCORE = e.score;

    return sf;
}

