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

#include "align.h"
#define MAX_TILE_SIZE 512

int INF =  std::numeric_limits<int>::max();
enum AlnOperands {ZERO_OP, DELETE_OP, INSERT_OP, MATCH_OP};
enum states {Z, D, I, M};

std::queue<int> AlignWithBT(char* ref_seq, long long int ref_len, char* query_seq, long long int query_len, int* sub_mat, int gap_open, int gap_extend, int query_pos, int ref_pos, bool reverse, bool first, int early_terminate) {

  // terminate if the R or Q length is greater than the tile size
  assert(ref_len < MAX_TILE_SIZE);
  assert(query_len < MAX_TILE_SIZE);
  
  int h_matrix_wr[MAX_TILE_SIZE + 1];
  int m_matrix_wr[MAX_TILE_SIZE + 1];
  int i_matrix_wr[MAX_TILE_SIZE + 1];
  int d_matrix_wr[MAX_TILE_SIZE + 1];

  int h_matrix_rd[MAX_TILE_SIZE + 1];
  int m_matrix_rd[MAX_TILE_SIZE + 1];
  int i_matrix_rd[MAX_TILE_SIZE + 1];
  int d_matrix_rd[MAX_TILE_SIZE + 1];

  int dir_matrix[MAX_TILE_SIZE+1][MAX_TILE_SIZE+1];


  for (int i = 0; i < query_len + 1; i++) {
    h_matrix_rd[i] = 0;
    m_matrix_rd[i] = 0;
    i_matrix_rd[i] = -INF;
    d_matrix_rd[i] = -INF;
   
    h_matrix_wr[i] = 0;
    m_matrix_wr[i] = 0;
    i_matrix_wr[i] = -INF;
    d_matrix_wr[i] = -INF;
  }
  
  for (int i = 0; i < ref_len + 1; i++) {
      dir_matrix[i][0] = ZERO_OP;
  }

  for (int j = 0; j < query_len + 1; j++) {
      dir_matrix[0][j] = ZERO_OP;
  }

  
  int max_score = 0; 
  int pos_score = 0; 
  int max_i = 0; 
  int max_j = 0; 
  
  for (int i = 1; i < ref_len + 1; i++) {
      for (int k = 1; k < MAX_TILE_SIZE + 1; k++) {
          m_matrix_rd[k] = m_matrix_wr[k];
          h_matrix_rd[k] = h_matrix_wr[k];
          i_matrix_rd[k] = i_matrix_wr[k];
          d_matrix_rd[k] = d_matrix_wr[k];
      }

      //j - row number; i - column number
      for (int j = 1; j < query_len + 1; j++) {
          int ref_nt = (reverse) ? NtChar2Int(ref_seq[ref_len-i]) : NtChar2Int(ref_seq[i-1]);
          int query_nt = (reverse) ? NtChar2Int(query_seq[query_len-j]) : NtChar2Int(query_seq[j-1]);

          int match;
          //case of unknown nucleotide in either reference or query
          match = sub_mat[query_nt*5 + ref_nt];

          //columnwise calculations
          m_matrix_wr[j] = h_matrix_rd[j-1] + match;

          int ins_open   = m_matrix_rd[j] + gap_open;
          int ins_extend = i_matrix_rd[j] + gap_extend;
          int del_open   = m_matrix_wr[j-1] + gap_open;
          int del_extend = d_matrix_wr[j-1] + gap_extend;

          i_matrix_wr[j] = (ins_open > ins_extend) ? ins_open : ins_extend;

          d_matrix_wr[j] = (del_open > del_extend) ? del_open : del_extend;

          int max1 = m_matrix_wr[j] > i_matrix_wr[j] ? m_matrix_wr[j] : i_matrix_wr[j];
          int max2 = d_matrix_wr[j] > 0 ? d_matrix_wr[j] : 0;
          h_matrix_wr[j] = max1 > max2 ? max1 : max2; 

          (dir_matrix)[i][j] = ((m_matrix_wr[j] >= i_matrix_wr[j]) ? ((m_matrix_wr[j] >= d_matrix_wr[j]) ? MATCH_OP : DELETE_OP) : ((i_matrix_wr[j] >= d_matrix_wr[j]) ? INSERT_OP : DELETE_OP));

          if ((m_matrix_wr[j] <= 0) && (i_matrix_wr[j] <= 0) && (d_matrix_wr[j] <= 0)) {
              (dir_matrix)[i][j] = ZERO_OP;
          }

          (dir_matrix)[i][j] += (ins_open >= ins_extend) ? (2 << INSERT_OP) : 0;
          (dir_matrix)[i][j] += (del_open >= del_extend) ? (2 << DELETE_OP) : 0;

          if (h_matrix_wr[j] >= max_score) {
              max_score = h_matrix_wr[j];
              max_i = i;
              max_j = j;
          }

          if ((i == ref_pos) && (j == query_pos)) {
              pos_score = h_matrix_wr[j];
          }

      }
  }

  std::queue<int> BT_states;

  int i_curr=ref_pos, j_curr=query_pos;
  int i_steps = 0, j_steps = 0;

  int open = 0;
  if (first) {
      i_curr = max_i;
      j_curr = max_j;
      BT_states.push(max_score);
      BT_states.push(i_curr);
      BT_states.push(j_curr);
  }
  else {
      BT_states.push(pos_score);
  }

  int state = dir_matrix[i_curr][j_curr] % 4;

  while (state != Z) {
      if ((i_steps >= early_terminate) || (j_steps >= early_terminate)) { 
          break;
      }
      BT_states.push(state);
      if (state == M) {
          state = (dir_matrix[i_curr-1][j_curr-1] % 4);
          i_curr--;
          j_curr--;
          i_steps++;
          j_steps++;
      }
      else if (state == I) {
          state = (dir_matrix[i_curr][j_curr] & (2 << INSERT_OP)) ? M : I;
          i_curr--;
          i_steps++;
      }
      else if (state == D) {
          state = (dir_matrix[i_curr][j_curr] & (2 << DELETE_OP)) ? M : D;
          j_curr--;
          j_steps++;
      }
  }; 

  return BT_states;
}


