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

`timescale 1ns / 1ps

module tb_GACTTop();
  localparam test_cycles = 50000;
  parameter HalfClkPeriod = 5;
  parameter REF_FILENAME = "../test_data/ref.10.hex";
  parameter QUERY_FILENAME = "../test_data/query.10.hex";
  
  localparam ClkPeriod = 2*HalfClkPeriod;
  reg [31:0] cycle;

  reg           clk;
  reg           rst;
  wire          done;
  
  wire [8:0] num_query_bases;
  wire [8:0] num_ref_bases;
  wire [17:0] num_tb_steps;
  wire [8:0] query_max_pos;
  wire       ready;
  wire [8:0] ref_max_pos;
  wire [15:0] req_id_out;
  wire [9:0]  tile_score;
  wire [1:0] dir;
  wire  dir_valid;

  reg clear_done;
  reg [8:0] ref_len;
  reg [8:0] query_len;

  reg [9:0] max_tb_steps;

  reg [9:0]  score_threshold;
  reg [119:0] in_params;

  reg [9:0] params [0:11];
  reg [7:0] align_fields;
  reg [15:0] req_id_in;

  reg set_params;
  reg start;

  reg [63:0] ref_in;
  reg [63:0] query_in;
  reg ref_wr_en;
  reg query_wr_en;
  reg [5:0] ref_addr_in; 
  reg [5:0] query_addr_in; 

  // Generate Clock
  initial clk = 1;
  always #(HalfClkPeriod) clk = ~clk;

  GACTTop #(
      .PE_WIDTH(10),
      .BLOCK_WIDTH(3),
      .MAX_TILE_SIZE(512),
      .NUM_PE(4),
      .REF_FILENAME(REF_FILENAME),
      .QUERY_FILENAME(QUERY_FILENAME),
      .REQUEST_ID_WIDTH(16)
  ) dut (
      .clk(clk),
      .rst(rst),
      .align_fields (align_fields),

      .clear_done (clear_done),
      .in_params (in_params),
      .max_tb_steps (max_tb_steps),
      .query_addr_in (query_addr_in),
      .query_in (query_in),
      .query_len (query_len),
      .query_wr_en (query_wr_en),
      .ref_addr_in (ref_addr_in),
      .ref_in (ref_in),
      .ref_len (ref_len),
      .ref_wr_en (ref_wr_en),
      .req_id_in (req_id_in),
      .score_threshold (score_threshold),
      .set_params (set_params),
      .start (start),

      .dir (dir),
      .dir_valid (dir_valid),

      .done (done),
      .num_query_bases (num_query_bases),
      .num_ref_bases (num_ref_bases),
      .num_tb_steps (num_tb_steps),
      .query_max_pos (query_max_pos),
      .ready (ready),
      .ref_max_pos (ref_max_pos),
      .req_id_out (req_id_out),
      .tile_score (tile_score)

  );

  // Run simulation 
  initial begin
      $dumpfile("test.10.vcd");
      $dumpvars(0, dut); 
      cycle = 0;
      start = 0;
      set_params = 0;

      score_threshold = 0;
      clear_done = 1;
      
      ref_len = 320;
      query_len = 320;
      ref_wr_en = 0;
      query_wr_en = 0;

      max_tb_steps = 400;
      
      align_fields = (1<<5);

      params[11] = 1;
      params[10] = -1;
      params[9] = -1;
      params[8] = -1;

      params[7] = 1;
      params[6] = -1;
      params[5] = -1;

      params[4] = 1;
      params[3] = -1;

      params[2] = 1;

      params[1] = -1;
      params[0] = -1;

      in_params = {params[11], params[10], params[9], params[8], params[7], params[6], params[5], params[4], params[3], params[2], params[1], params[0]};            
      $display("---- Performing Reset ----");
      rst = 1; 
      
      #(1*ClkPeriod); 
      rst = 0;
      
      #(1*ClkPeriod); 
      set_params = 1;

      #(1*ClkPeriod); 
      set_params = 0;
      
      #(1*ClkPeriod); 
      $display("---- Starting Simulation ----");
      start = 1;

      #(1*ClkPeriod); 
      start = 0;
  end
  
  always @ (posedge clk) begin
      cycle <= cycle + 1;
      if (cycle > test_cycles) begin
          $display("---- Done ----");
          $finish();
      end
  end

  always @ (posedge clk) begin
      if (dir_valid) begin
          $display("%d", dir);
      end
      if (done) begin
          $display("ref pos: %d, query pos: %d, max score", ref_max_pos, query_max_pos, tile_score);
      end
  end

endmodule
