 /*
MIT License

Copyright (c) 2018 Yatish Turakhia, Sneha D. Goenka, Gill Bejerano and William Dally

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

module BTLogic #(
  parameter ADDR_WIDTH = 20,
  parameter REF_LEN_WIDTH = 12,
  parameter LOG_NUM_PE = 6
)
(
    input clk,
    input rst,
    input start,

    input [REF_LEN_WIDTH-1:0] ref_length,
    input [ADDR_WIDTH-1:0] max_score_mod_addr,
    input [ADDR_WIDTH-1:0] max_score_addr,
    input [LOG_NUM_PE-1:0] max_score_pe,
    input [1:0] max_score_pe_state,
    input [3:0] input_dir,
    input [3:0] input_dir_diag,

    output reg [ADDR_WIDTH-1:0] next_addr,
    output reg [LOG_NUM_PE-1:0] next_pe,
    output wire [ADDR_WIDTH-1:0] next_addr_diag,
    output wire [LOG_NUM_PE-1:0] next_pe_diag,
    output wire addr_valid,
    output wire [1:0] dir,
    output wire dir_valid,
    output reg [REF_LEN_WIDTH-1:0] H_offset,
    input [REF_LEN_WIDTH-1:0] max_H_offset,
    output reg [REF_LEN_WIDTH-1:0] V_offset,
    input [REF_LEN_WIDTH-1:0] max_V_offset,
    output reg [ADDR_WIDTH+LOG_NUM_PE-1:0] num_tb_steps,
    output done
);

  localparam MAX_PE = (2**LOG_NUM_PE) - 1;

  reg [ADDR_WIDTH-1:0] mod_count;
  reg [1:0] next_pe_state;

  reg [ADDR_WIDTH-1:0] mod_next_addr;

  reg [REF_LEN_WIDTH-1:0] reg_ref_length;

  reg [2:0] state;
  reg [2:0] next_state;
  
  localparam WAIT=0, BLOCK1=1, BLOCK2=2, CALC=3, DONE=4;
  localparam ZERO=0, M=3, V=1, H=2;

  wire [LOG_NUM_PE-1:0] next_pe_decr_mod;
  wire [ADDR_WIDTH-1:0] next_addr_mod;
  
  wire [REF_LEN_WIDTH-1:0] next_V_offset;
  wire [REF_LEN_WIDTH-1:0] next_H_offset;

  assign done = (state == DONE);
  assign dir_valid = (state == CALC) && (dir != 0);
  assign addr_valid = (next_state == CALC);
  assign dir = next_pe_state;

  assign next_pe_decr_mod = (next_pe == 0) ? MAX_PE : (next_pe - 1);
  assign next_addr_mod = (next_pe == 0) ? (next_addr - reg_ref_length - 1) : (next_addr - 1);
  
  assign next_addr_diag = next_addr_mod;
  assign next_pe_diag = next_pe_decr_mod;
  
  assign next_V_offset = ((next_pe_state == M) || (next_pe_state == V)) ? (V_offset + 1) : V_offset;
  assign next_H_offset = ((next_pe_state == M) || (next_pe_state == H)) ? (H_offset + 1) : H_offset;

  always @(posedge clk) begin
      mod_next_addr <= mod_count; //??
  end

  always @(posedge clk) begin
      if (rst) begin
          state <= WAIT;
      end
      else begin
          state <= next_state;
          case (state)
              WAIT: begin
                  mod_count <= max_score_mod_addr;
                  next_addr <= max_score_addr;
                  next_pe <= max_score_pe; 
                  next_pe_state <= max_score_pe_state;
                  reg_ref_length <= ref_length;
                  H_offset <= 0;
                  V_offset <= 0;
                  num_tb_steps <= 0;
              end
              BLOCK1: begin
              end
              BLOCK2: begin
              end
              CALC: begin
                  H_offset <= next_H_offset;
                  V_offset <= next_V_offset;
                  num_tb_steps <= num_tb_steps + (dir != 0);
                  if (next_pe_state == M) begin
                      //terminating condition
                      if (((next_pe == 0) && (next_addr <= reg_ref_length)) || (mod_next_addr == 0)) begin
                          next_pe_state <= ZERO;
                      end
                      //otherwise
                      else begin
                          mod_count <= mod_count - 1;
                          if (next_pe == 0) begin
                              next_addr <= next_addr - reg_ref_length - 1;
                              next_pe <= MAX_PE;
                          end
                          else begin
                              next_addr <= next_addr - 1;
                              next_pe <= next_pe - 1;
                          end
                          next_pe_state <= input_dir_diag[1:0];
                      end
                 end
                 else if (next_pe_state == V) begin
                     if (next_pe == 0) begin
                         next_pe <= MAX_PE;
                     end
                     else begin
                         next_pe <= next_pe - 1;
                     end
                     if ((next_pe == 0) && (next_addr <= reg_ref_length)) begin
                         next_pe_state <= ZERO;
                     end
                     else begin
                         if (next_pe == 0) begin
                             next_addr <= next_addr - reg_ref_length;
                         end
                         if (input_dir[2] == 0) begin
                             next_pe_state <= V;
                         end
                         else begin 
                             next_pe_state <= M;
                         end
                     end
                 end
                 else if (next_pe_state == H) begin
                     next_addr <= next_addr - 1;
                     mod_count <= mod_count - 1;
                     if (mod_count == 0) begin 
                         next_pe_state <= ZERO;
                     end
                     else begin
                         if (input_dir[3] == 0) begin
                             next_pe_state <= H;
                         end
                         else begin
                             next_pe_state <= M;
                         end
                     end
                 end
             end
             DONE: begin
             end
         endcase
     end
  end

  always @(*) 
  begin
      next_state = state;
      case (state)
          WAIT:
              if (start)
                  next_state = BLOCK1;
          BLOCK1:
              next_state = BLOCK2;
          BLOCK2:
              next_state = CALC;
          CALC:
              if ((next_pe_state == ZERO) || (next_H_offset == max_H_offset) || (next_V_offset == max_V_offset))
                  next_state = DONE;
              else
                  next_state = BLOCK1;
          DONE:
              next_state = WAIT;
      endcase
  end

endmodule

