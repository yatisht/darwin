/*
MIT License

Copyright (c) 2018 Yatish Turakhia, Gill Bejerano and William Dally

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "gact.h"

enum states {Z, D, I, M};


Alignment GACT (char* ref_str, char* query_str, std::string ref_name, std::string query_name, int* sub_mat, int gap_open, int gap_extend, int tile_size, int tile_overlap, int ref_pos, int query_pos, uint32_t ref_length, uint32_t query_length, char strand, int first_tile_score_threshold) {
    std::queue<int> BT_states;
   
    std::string aligned_ref_str = "";
    std::string aligned_query_str = "";

    Alignment alignment;
    //output of the function
    alignment.ref_name = ref_name;
    alignment.query_name = query_name;
    alignment.aligned_ref_str = "";
    alignment.aligned_query_str = "";
    alignment.ref_start = ref_pos;
    alignment.query_start = query_pos;
    alignment.aligned_ref_len = 0;
    alignment.aligned_query_len = 0;
    alignment.ref_len = ref_length;
    alignment.query_len = query_length;
    alignment.strand = strand;
    alignment.score = 0;
    alignment.flag = 0;
    
    int ref_tile_length = tile_size;
    int query_tile_length = tile_size;
    
//    //length of the complete sequences
//    int ref_length = ref_str.length();
//    int query_length = query_str.length();
    
    // beginning for the left or the right extension
    //original starting point for the whole query
    int rev_ref_pos = ref_pos;
    int rev_query_pos = query_pos;
   
    int max_ref_pos = 0;
    int max_query_pos = 0;
    int i = 0;
    int j = 0;
    
    int first_tile_score = 0;
    bool first_tile = true;
   
    // not for the first tile
    while ((ref_pos > 0) && (query_pos > 0) && (((i > 0) && (j > 0)) || first_tile)) {
        //change the tile length if elements less than that of the tile size
        ref_tile_length = (ref_pos > tile_size) ? tile_size : ref_pos;
        query_tile_length = (query_pos > tile_size) ? tile_size : query_pos;

//        std::cout << std::string(ref_str+(ref_pos-ref_tile_length), ref_tile_length) << "\n" << std::string(query_str+(query_pos-query_tile_length), query_tile_length) << "\n";
        BT_states = AlignWithBT (ref_str+(ref_pos-ref_tile_length), ref_tile_length, query_str+(query_pos-query_tile_length), query_tile_length, sub_mat, gap_open, gap_extend, query_tile_length, ref_tile_length, false, first_tile, (tile_size - tile_overlap));
        i = 0;
        j = 0;
        int tile_score = BT_states.front();
        BT_states.pop();
        
        if (first_tile) {
            ref_pos = ref_pos - ref_tile_length + BT_states.front();
            max_ref_pos = BT_states.front();
            BT_states.pop();
            query_pos = query_pos - query_tile_length + BT_states.front();
            max_query_pos = BT_states.front();
            BT_states.pop();
            rev_ref_pos = ref_pos;
            rev_query_pos = query_pos;
            first_tile_score = tile_score;
//            if (tile_score < first_tile_score_threshold) {
//                break;
//            }
        }

        int num_tb = BT_states.size();
        char* ref_buf = (char*) malloc(num_tb);
        char* query_buf = (char*) malloc(num_tb);
        int ref_buf_curr = num_tb-1;
        int query_buf_curr = num_tb-1;

        while (!BT_states.empty()) {
            first_tile = false;
            int state = BT_states.front();
            BT_states.pop();
            if (state == M) {
//                aligned_ref_str.insert(0, 1, ref_str[ref_pos - j - 1]);
//                aligned_query_str.insert(0, 1, query_str[query_pos - i - 1]);
                ref_buf[ref_buf_curr--] = ref_str[ref_pos - j - 1];
                query_buf[query_buf_curr--] = query_str[query_pos - i - 1];
                i += 1;
                j += 1;
            }
            if (state == I) {
//                aligned_ref_str.insert(0, 1, ref_str[ref_pos - j - 1]);
//                aligned_query_str.insert(0, 1, '-');
                ref_buf[ref_buf_curr--] = ref_str[ref_pos - j - 1];
                query_buf[query_buf_curr--] = '-';
                j += 1;
            }
            if (state == D) {
//                aligned_ref_str.insert(0, 1, '-');
//                aligned_query_str.insert(0, 1, query_str[query_pos - i - 1]);
                ref_buf[ref_buf_curr--] = '-';
                query_buf[query_buf_curr--] = query_str[query_pos - i - 1];
                i += 1;
            }
        }
        if (num_tb > 0) {
            aligned_ref_str = std::string(ref_buf, num_tb) + aligned_ref_str;
            aligned_query_str = std::string(query_buf, num_tb) + aligned_query_str;
        }
        free(ref_buf);
        free(query_buf);
        ref_pos -= (j);
        query_pos -= (i);
        alignment.aligned_ref_len += j;
        alignment.aligned_query_len += i;
    }
    
    alignment.ref_start = ref_pos;
    alignment.query_start = query_pos;

    ref_pos = rev_ref_pos;
    query_pos = rev_query_pos;
    
    i =  tile_size;
    j = tile_size;
    
    //starts with the first tile
    while ((ref_pos < ref_length) && (query_pos < query_length) && (((i > 0) && (j > 0)) || first_tile)) {
        ref_tile_length = (ref_pos + tile_size < ref_length) ? tile_size : ref_length - ref_pos;
        query_tile_length = (query_pos + tile_size < query_length) ? tile_size : query_length - query_pos;
        BT_states = AlignWithBT (ref_str+ref_pos, ref_tile_length, query_str+query_pos, query_tile_length, sub_mat, gap_open, gap_extend, query_tile_length, ref_tile_length, true, first_tile, (tile_size - tile_overlap));
        i = 0;
        j = 0;
        int tile_score = BT_states.front();
        BT_states.pop();
        if (first_tile) {
            ref_pos = ref_pos + ref_tile_length - BT_states.front();
            max_ref_pos = BT_states.front();
            BT_states.pop();
            query_pos = query_pos + query_tile_length - BT_states.front();
            max_query_pos = BT_states.front();
            BT_states.pop();
            first_tile_score = tile_score;
//            if (tile_score < first_tile_score_threshold) {
//                break;
//            }
        }

        int num_tb = BT_states.size();
        char* ref_buf = (char*) malloc(num_tb);
        char* query_buf = (char*) malloc(num_tb);
        int ref_buf_curr = 0;
        int query_buf_curr = 0;
        
        while (!BT_states.empty()) {
            first_tile = false;
            int state = BT_states.front();
            BT_states.pop();
            if (state == M) {
//                aligned_ref_str += ref_str[ref_pos + j];
//                aligned_query_str += (query_str[query_pos + i]);
                ref_buf[ref_buf_curr++] = ref_str[ref_pos + j];
                query_buf[query_buf_curr++] = query_str[query_pos + i];
                i += 1;
                j += 1;
            }
            if (state == I) {
//                aligned_ref_str += ref_str[ref_pos + j];
//                aligned_query_str += '-';
                ref_buf[ref_buf_curr++] = ref_str[ref_pos + j];
                query_buf[query_buf_curr++] = '-';
                j += 1;
            }
            if (state == D) {
//                aligned_ref_str += '-';
//                aligned_query_str += query_str[query_pos + i];
                ref_buf[ref_buf_curr++] = '-';
                query_buf[query_buf_curr++] = query_str[query_pos + i];
                i += 1;
            }
        }
        if (num_tb > 0) {
            aligned_ref_str = aligned_ref_str + std::string(ref_buf, num_tb);
            aligned_query_str = aligned_query_str + std::string(query_buf, num_tb);
        }
        free(ref_buf);
        free(query_buf);
        ref_pos += (j);
        query_pos += (i);
        alignment.aligned_ref_len += j;
        alignment.aligned_query_len += i;
    }

    int total_score = 0;
    bool open = true;
    for (uint32_t j = 0; j < aligned_ref_str.length(); j++) {
        char ref_nt = aligned_ref_str[j];
        char query_nt = aligned_query_str[j];
        if (ref_nt == '-' || query_nt == '-') {
            total_score += (open) ? gap_open : gap_extend;
            open = false;
        }
        else {
            total_score += sub_mat[5*NtChar2Int(query_nt) + NtChar2Int(ref_nt)];
            open = true;
        }
    }
    alignment.aligned_ref_str = aligned_ref_str;
    alignment.aligned_query_str = aligned_query_str;
    alignment.score = total_score;
//    std::cout << aligned_ref_str << std::endl << aligned_query_str << std::endl;
//    std::cout << max_ref_pos-1 << " " << max_query_pos-1 << std::endl;
//    std::cout << "Max score: " <<  first_tile_score << std::endl;
//    std::cout << "Total score: " << total_score << std::endl;
//    std::cout  << std::endl;
    return alignment;
}

