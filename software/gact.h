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

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cstdlib>
#include <queue>
#include "align.h"

struct Alignment {
    std::string ref_name;
    std::string query_name;
    std::string aligned_ref_str;
    std::string aligned_query_str;
    uint32_t ref_start;
    uint32_t query_start;
    uint32_t aligned_ref_len;
    uint32_t aligned_query_len;
    uint32_t ref_len;
    uint32_t query_len;
    int score;
    int flag;
    char strand;
};

Alignment GACT (char* ref_str, char* query_str, std::string ref_name, std::string query_name, int* sub_mat, int gap_open, int gap_extend, int tile_size, int tile_overlap, int ref_pos, int query_pos, uint32_t ref_length, uint32_t query_length, char strand, int first_tile_score_threshold);
