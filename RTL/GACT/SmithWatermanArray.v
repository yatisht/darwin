
/*
MIT License

Copyright (c) 2018 Yatish Turakhia, Sneha D. Goenka, Albert Ng, Gill Bejerano and William Dally

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

/* ref_length should be greater than NUM_PE
* and NUM_PE should be power of 2.
*/
module SmithWatermanArray #(
    parameter NUM_PE = 32,
    parameter LOG_NUM_PE = 5,
    parameter REF_LEN_WIDTH = 10,
    parameter QUERY_LEN_WIDTH = 10,
    parameter PE_WIDTH = 8,
    parameter PARAM_ADDR_WIDTH = 10,
    
    parameter REF_BLOCK_SIZE_WIDTH = 7, // REF_LEN_WIDTH;
    parameter QUERY_BLOCK_SIZE_WIDTH = 7 //QUERY_LEN_WIDTH;
)(
    input  clk,         
    input  rst,        
    input start,
    input [12*PE_WIDTH-1:0] in_param,

    input reverse_ref_in,
    input reverse_query_in,

    input complement_ref_in,
    input complement_query_in,

    input do_traceback_in,
    input [REF_LEN_WIDTH-1:0] ref_length,
    input [QUERY_LEN_WIDTH-1:0] query_length,

    input [REF_LEN_WIDTH-1:0] ref_bram_rd_start_addr, 
    output reg [REF_LEN_WIDTH-1:0] ref_bram_rd_addr, 
    input [7:0] ref_bram_data_in,
    input [PARAM_ADDR_WIDTH-1:0] query_bram_rd_start_addr,
    output reg [PARAM_ADDR_WIDTH-1:0] query_bram_rd_addr,
    input [7:0] query_bram_data_in,

    input [PE_WIDTH-1:0] max_score_threshold,
    input start_last,

    output reg [PE_WIDTH-1:0] max_score,
    output reg [REF_LEN_WIDTH-1:0] H_offset,
    input [REF_BLOCK_SIZE_WIDTH-1:0] max_H_offset,
    output reg [QUERY_LEN_WIDTH-1:0] V_offset,
    input [QUERY_BLOCK_SIZE_WIDTH-1:0] max_V_offset,

    output reg [REF_LEN_WIDTH-1:0] ref_max_score_pos,
    output reg [QUERY_LEN_WIDTH-1:0] query_max_score_pos,

    output reg [(REF_BLOCK_SIZE_WIDTH + (QUERY_BLOCK_SIZE_WIDTH - LOG_NUM_PE))+LOG_NUM_PE-1:0] num_tb_steps,

    output wire [1:0] dir,
    output wire dir_valid,

    output done);

    localparam PARAM_WIDTH = 7 * PE_WIDTH;
    
    localparam BT_BRAM_ADDR_WIDTH = REF_BLOCK_SIZE_WIDTH + (QUERY_BLOCK_SIZE_WIDTH - LOG_NUM_PE);
    localparam V_FIFO_DEPTH_WIDTH = REF_BLOCK_SIZE_WIDTH;
    localparam M_FIFO_DEPTH_WIDTH = REF_BLOCK_SIZE_WIDTH;
    localparam F_FIFO_DEPTH_WIDTH = REF_BLOCK_SIZE_WIDTH;

    reg first_query_block;

    reg ref_fifo_rd_en;

    reg [NUM_PE-1:0] set_param;
    reg [NUM_PE-1:0] last;
    wire [PARAM_WIDTH-1:0] param;
    wire [PARAM_WIDTH-1:0] param_out;

    reg [REF_LEN_WIDTH-1:0] curr_ref_len;
    reg [QUERY_LEN_WIDTH-1:0] curr_query_len;
    
    reg reverse_query, reverse_ref;
    reg complement_query, complement_ref;
    reg do_traceback;

    reg rst_pe[NUM_PE-1:0];
    reg [3:0] state;
    reg [3:0] next_state;

    wire [3:0] ref_nt; 
    wire [3:0] query_nt; 

    wire [4*NUM_PE-1:0] select_dir_data_in;
  
    wire [PE_WIDTH-1:0] sub_A_in;
    wire [PE_WIDTH-1:0] sub_C_in;
    wire [PE_WIDTH-1:0] sub_G_in;
    wire [PE_WIDTH-1:0] sub_T_in;
    wire [PE_WIDTH-1:0] sub_N_in;
    wire [PE_WIDTH-1:0] gap_extend_in;
    wire [PE_WIDTH-1:0] gap_open_in;

    wire [PE_WIDTH-1:0] F_in[0:NUM_PE-1];
    wire [2:0] T_in[0:NUM_PE-1];
    wire [PE_WIDTH-1:0] V_in[0:NUM_PE-1];
    wire [PE_WIDTH-1:0] M_in[0:NUM_PE-1];
    wire [PE_WIDTH-1:0] E_in[0:NUM_PE-1];
    wire [PE_WIDTH-1:0] F_out[0:NUM_PE-1];
    wire [2:0] T_out[0:NUM_PE-1];
    wire [PE_WIDTH-1:0] V_out[0:NUM_PE-1];
    wire [PE_WIDTH-1:0] M_out[0:NUM_PE-1];
    wire [PE_WIDTH-1:0] E_out[0:NUM_PE-1];

    wire [PE_WIDTH-1:0] init_V;
    wire [PE_WIDTH-1:0] init_E;

    wire [3:0] pe_dir [0:NUM_PE-1];
    wire [BT_BRAM_ADDR_WIDTH-1:0] pe_dir_addr [0:NUM_PE-1];
    wire pe_dir_valid [0:NUM_PE-1];

    wire init_in [0:NUM_PE-1];
    wire init_out [0:NUM_PE-1];

    wire [BT_BRAM_ADDR_WIDTH-1:0] max_ref_pos_in [0:NUM_PE-1];
    wire [BT_BRAM_ADDR_WIDTH-1:0] max_ref_pos_out [0:NUM_PE-1];
    wire [BT_BRAM_ADDR_WIDTH-1:0] max_ref_mod_in [0:NUM_PE-1];
    wire [BT_BRAM_ADDR_WIDTH-1:0] max_ref_mod_out [0:NUM_PE-1];
    wire [LOG_NUM_PE-1:0] max_query_pos_in [0:NUM_PE-1];
    wire [LOG_NUM_PE-1:0] max_query_pos_out [0:NUM_PE-1];
    wire [QUERY_BLOCK_SIZE_WIDTH-1:0] max_query_mod_in [0:NUM_PE-1];
    wire [QUERY_BLOCK_SIZE_WIDTH-1:0] max_query_mod_out [0:NUM_PE-1];
    wire [1:0] max_pe_state_in [0:NUM_PE-1];
    wire [1:0] max_pe_state_out [0:NUM_PE-1];

    wire compute_max_in [0:NUM_PE-1];
    wire compute_max_out [0:NUM_PE-1];
    wire last_in [0:NUM_PE-1];
    wire last_out [0:NUM_PE-1];
    wire stall;

    wire [BT_BRAM_ADDR_WIDTH-1:0] bt_bram_addr[0:NUM_PE-1];
    wire [3:0] bt_bram_data_out[0:NUM_PE-1];

    wire V_intern_fifo_wr_en;
    wire V_intern_fifo_rd_en;
    wire [PE_WIDTH-1:0] V_intern_fifo_wr_data_in;
    wire [PE_WIDTH-1:0] V_intern_fifo_rd_data_out;
    wire V_intern_fifo_full;
    wire V_intern_fifo_empty;
  
    wire M_intern_fifo_wr_en;
    wire M_intern_fifo_rd_en;
    wire [PE_WIDTH-1:0] M_intern_fifo_wr_data_in;
    wire [PE_WIDTH-1:0] M_intern_fifo_rd_data_out;
    wire M_intern_fifo_full;
    wire M_intern_fifo_empty;
  
    wire F_intern_fifo_wr_en;
    wire F_intern_fifo_rd_en;
    wire [PE_WIDTH-1:0] F_intern_fifo_wr_data_in;
    wire [PE_WIDTH-1:0] F_intern_fifo_rd_data_out;
    wire F_intern_fifo_full;
    wire F_intern_fifo_empty;
  
    wire bt_logic_start;
    wire [BT_BRAM_ADDR_WIDTH-1:0] bt_logic_max_score_addr;
    wire [BT_BRAM_ADDR_WIDTH-1:0] bt_logic_max_score_mod_addr;
    wire [LOG_NUM_PE-1:0] bt_logic_max_score_pe;
    wire [1:0] bt_logic_max_score_pe_state;
    wire [3:0] bt_logic_input_dir;
    wire [3:0] bt_logic_input_dir_diag;
    wire [BT_BRAM_ADDR_WIDTH-1:0] bt_logic_next_addr;
    wire [BT_BRAM_ADDR_WIDTH-1:0] bt_logic_next_addr_diag;
    wire [LOG_NUM_PE-1:0] bt_logic_next_pe;
    wire [LOG_NUM_PE-1:0] bt_logic_next_pe_diag;
    wire bt_logic_addr_valid;
    wire [1:0] bt_logic_dir_out;
    wire [REF_BLOCK_SIZE_WIDTH-1:0] bt_logic_H_offset;
    wire [REF_BLOCK_SIZE_WIDTH-1:0] bt_logic_V_offset;
    wire bt_logic_dir_valid;
    wire bt_logic_done;
    reg [REF_BLOCK_SIZE_WIDTH-1:0] bt_ref_length;
    wire[BT_BRAM_ADDR_WIDTH+LOG_NUM_PE-1:0] bt_logic_num_tb_steps;

    reg [PE_WIDTH-1:0] F_in_0;
    reg [PE_WIDTH-1:0] V_in_0;
    reg [PE_WIDTH-1:0] M_in_0;
    reg [3:0] T_in_0;
    reg init_in_0;
    reg delayed_init_in_0;
    reg compute_max_in_0;
    reg [PE_WIDTH-1:0] gap_open;
    reg [PE_WIDTH-1:0] gap_extend;

    Ascii2Nt ref_ascii2nt
    (
        .ascii(ref_bram_data_in),
        .complement(complement_ref),
        .nt(ref_nt)
    );

    Ascii2Nt query_ascii2nt (
        .ascii(query_bram_data_in),
        .complement(complement_query),
        .nt(query_nt)
    );

    Nt2Param #(
        .PE_WIDTH(PE_WIDTH)
    ) query_nt2param (
        .nt(query_nt),
        .in_param(in_param),
        .out_param(param_out)
    );

    genvar k;

    localparam WAIT=0, READ_PARAM_FIFO=1, SET_PARAM=2, STREAM_REF_START=3, STREAM_REF=4, STREAM_REF_STOP=5, STREAM_REF_DONE=6, COMPUTE_MAX_START=7, COMPUTE_MAX_WAIT=8, BT_START=9, BT_WAIT=10, DONE=11;

    assign dir = bt_logic_dir_out; 
    assign dir_valid = bt_logic_dir_valid;

    assign param = (curr_query_len > query_length) ? 0 : param_out;

    assign done = (state == DONE);

  always @(posedge clk) begin
      if (rst) begin
          set_param <= 0;
          last <= 0;
          first_query_block <= 1;
          ref_fifo_rd_en <= 0;
          curr_query_len <= 0;
          H_offset <= 0;
          V_offset <= 0;
          init_in_0 <= 0;              
          state <= WAIT;
          gap_open <= in_param[2*PE_WIDTH-1-:PE_WIDTH];
          gap_extend <= in_param[PE_WIDTH-1:0];
      end
      else begin
          state <= next_state;
          delayed_init_in_0 <= init_in_0;

          if (compute_max_out[NUM_PE-1] == 1'b1) begin
              max_score <= V_out[NUM_PE-1];
          end
          if (state == WAIT) begin
              bt_ref_length <= ref_length;
              reverse_query <= reverse_query_in;
              reverse_ref <= reverse_ref_in;
              complement_query <= complement_query_in;
              complement_ref <= complement_ref_in;
              do_traceback <= do_traceback_in;
          end
          if (state == READ_PARAM_FIFO) begin
              if (reverse_query) begin
                  query_bram_rd_addr <= query_bram_rd_start_addr + query_length - 1;
              end
              else begin
                  query_bram_rd_addr <= query_bram_rd_start_addr;
              end
          end
          if (state == SET_PARAM) begin
              if (reverse_query) begin
                  query_bram_rd_addr <= query_bram_rd_addr - 1;
              end
              else begin
                  query_bram_rd_addr <= query_bram_rd_addr + 1;
              end
              curr_query_len <= curr_query_len + 1;
              if (set_param == 0) begin
                  set_param <= 1;
                  if ((start_last == 1'b1) && (curr_query_len + 1 == query_length)) begin
                      last <= 1;
                  end
                  else begin
                      last <= 0;
                  end
              end
              else begin
                  set_param <= (set_param << 1);
                  if ((start_last == 1'b1) && (curr_query_len + 1 == query_length)) begin
                      last <= (set_param << 1);
                  end
                  else begin
                      last <= 0;
                  end
              end
          end
          if (state == STREAM_REF_START) begin
              set_param <= 0;
              last <= 0;
              curr_ref_len <= 1;
              if (reverse_ref) begin
                  ref_bram_rd_addr <= ref_bram_rd_start_addr + ref_length - 1;
              end
              else begin
                  ref_bram_rd_addr <= ref_bram_rd_start_addr;
              end
              ref_fifo_rd_en <= 1;
          end
          if (state == STREAM_REF) begin
              init_in_0 <= 1;              
              curr_ref_len <= curr_ref_len + 1;
              if (reverse_ref) begin
                  ref_bram_rd_addr <= ref_bram_rd_addr - 1;
              end
              else begin
                  ref_bram_rd_addr <= ref_bram_rd_addr + 1;
              end
              if (curr_ref_len >= ref_length) begin
                  ref_fifo_rd_en <= 0;
              end
          end
          if (state == STREAM_REF_STOP) begin
              init_in_0 <= 0;              
              first_query_block <= 0;
          end
          if (state == BT_WAIT) begin
              if (bt_logic_done) begin
                  H_offset <= bt_logic_H_offset;
                  V_offset <= bt_logic_V_offset;
                  num_tb_steps <= bt_logic_num_tb_steps;
              end
          end
      end
  end

  assign stall = 0;

  genvar j;
      assign {sub_A_in, sub_C_in, sub_G_in, sub_T_in, sub_N_in, gap_open_in, gap_extend_in} = param;
      assign init_V = 0;
      assign init_E = ((2'b11) << (PE_WIDTH-2)); // hack

  assign compute_max_in[0] = compute_max_in_0;
  assign last_in[0] = 0;
  assign max_ref_pos_in[0] = 0;
  assign max_ref_mod_in[0] = 0;
  assign max_query_pos_in[0] = 0;
  assign max_query_mod_in[0] = 0;
  assign max_pe_state_in[0] = 0;
  assign F_in[0] = F_in_0;
  assign V_in[0] = V_in_0;
  assign M_in[0] = V_in_0;
  assign T_in[0] = T_in_0; 
  assign init_in[0] = delayed_init_in_0; 

  always @(posedge clk) begin
      if (rst) begin
          F_in_0 <= 0;
          V_in_0 <= 0;
          M_in_0 <= 0;
          T_in_0 <= 0;
          ref_max_score_pos <=0;
          query_max_score_pos <= 0;
          compute_max_in_0 <= 0;
      end
      else begin
          F_in_0 <= (first_query_block) ? gap_open_in : F_intern_fifo_rd_data_out; // ? 
          V_in_0 <= (first_query_block || (state == COMPUTE_MAX_START)) ? 0 : V_intern_fifo_rd_data_out; // ? 
          M_in_0 <= (first_query_block || (state == COMPUTE_MAX_START)) ? 0 : M_intern_fifo_rd_data_out; // ? 
          T_in_0 <= ref_nt; 
          compute_max_in_0 <= (state == COMPUTE_MAX_START);
          if (((next_state == BT_START) && do_traceback) || ((next_state == DONE) && ~do_traceback)) begin
              ref_max_score_pos <= ref_bram_rd_start_addr + max_ref_mod_out[NUM_PE-1];
              query_max_score_pos <= query_bram_rd_start_addr + (max_query_mod_out[NUM_PE-1] << LOG_NUM_PE) + max_query_pos_out[NUM_PE-1];
          end
      end
  end

  generate
  for (j = 0; j < NUM_PE; j = j+1)
  begin: rst_pe_gen
      always @(posedge clk) begin
          if (rst) begin
              rst_pe[j] <= 1;
          end
          else begin
              rst_pe[j] <= 0;
          end
      end
  end
  endgenerate

  generate
  for (j = 1; j < NUM_PE; j=j+1) 
  begin:systolic_connections
      assign F_in[j] = F_out[j-1]; 
      assign T_in[j] = T_out[j-1]; 
      assign M_in[j] = M_out[j-1]; 
      assign V_in[j] = V_out[j-1]; 
      assign compute_max_in[j] = compute_max_out[j-1];
      assign last_in[j] = last_out[j-1];
      assign max_ref_pos_in[j] = max_ref_pos_out[j-1];
      assign max_ref_mod_in[j] = max_ref_mod_out[j-1];
      assign max_query_pos_in[j] = max_query_pos_out[j-1]; 
      assign max_query_mod_in[j] = max_query_mod_out[j-1]; 
      assign max_pe_state_in[j] = max_pe_state_out[j-1]; 
      assign init_in[j] = init_out[j-1]; 
  end
  endgenerate
  
  generate
  for (k = 0; k < NUM_PE; k = k + 1) 
  begin:bram_gen
      BRAM #(
          .ADDR_WIDTH(BT_BRAM_ADDR_WIDTH),
          .DATA_WIDTH(4)
      ) bt_bram (
      .clk(clk),
      .addr(bt_bram_addr[k]),
      .write_en(pe_dir_valid[k]),
      .data_in(pe_dir[k]),
      .data_out(bt_bram_data_out[k]));
  end
  endgenerate


  //generating the PEs
  generate
  for (k = 0; k < NUM_PE; k = k + 1) 
  begin:swpe_gen 
  SmithWatermanPE #(
          .WIDTH(PE_WIDTH),
          .REF_LEN_WIDTH(BT_BRAM_ADDR_WIDTH),
          .QUERY_LEN_WIDTH(QUERY_BLOCK_SIZE_WIDTH),
          .LOG_NUM_PE(LOG_NUM_PE),
          .PE_ID(k)
      ) swpe (
          .clk (clk),
          .rst (rst_pe[k]),
          .set_param (set_param[k]),
          .last (last[k]),

          .F_in (F_in[k]),
          .T_in (T_in[k]),
          .M_in (M_in[k]),
          .V_in (V_in[k]),

          .gap_extend_in (gap_extend_in),
          .gap_open_in (gap_open_in),

          .init_E (init_E),
          .init_V (init_V),
          .init_in (init_in[k]),

          .sub_A_in (sub_A_in),
          .sub_C_in (sub_C_in),
          .sub_G_in (sub_G_in),
          .sub_N_in (sub_N_in),
          .sub_T_in (sub_T_in),

          .compute_max_in(compute_max_in[k]),
          .compute_max_out(compute_max_out[k]),

          .last_in (last_in[k]),
          .last_out (last_out[k]),
          
          .max_ref_pos_in (max_ref_pos_in[k]),
          .max_ref_mod_in (max_ref_mod_in[k]),
          .max_ref_pos_out (max_ref_pos_out[k]),
          .max_ref_mod_out (max_ref_mod_out[k]),
          .max_query_pos_in (max_query_pos_in[k]),
          .max_query_mod_in (max_query_mod_in[k]),
          .max_query_pos_out (max_query_pos_out[k]),
          .max_query_mod_out (max_query_mod_out[k]),
          .max_pe_state_in (max_pe_state_in[k]),
          .max_pe_state_out (max_pe_state_out[k]),

          .E_out (E_out[k]),
          .F_out (F_out[k]),
          .T_out (T_out[k]),
          .M_out (M_out[k]),
          .V_out (V_out[k]),
          .dir (pe_dir[k]),
          .dir_addr (pe_dir_addr[k]),
          .dir_valid (pe_dir_valid[k]),
          .init_out (init_out[k]));
  end
  endgenerate

  assign V_intern_fifo_wr_en = init_out[NUM_PE-1];
  assign V_intern_fifo_wr_data_in = V_out[NUM_PE-1];
  assign V_intern_fifo_rd_en = ((ref_fifo_rd_en == 1'b1) && (first_query_block == 1'b0));

  FIFOWithCount #(
      .DATA_WIDTH(PE_WIDTH),
      .DEPTH_WIDTH(V_FIFO_DEPTH_WIDTH)
  ) V_intern_fifo
  (
      .clk(clk),
      .rst(rst),

      .wr_en(V_intern_fifo_wr_en),
      .rd_en(V_intern_fifo_rd_en),
      .wr_data_in(V_intern_fifo_wr_data_in),

      .rd_data_out(V_intern_fifo_rd_data_out),
      .count(),
      .full(V_intern_fifo_full),
      .empty(V_intern_fifo_empty)
  );

  assign M_intern_fifo_wr_en = init_out[NUM_PE-1];
  assign M_intern_fifo_wr_data_in = M_out[NUM_PE-1];
  assign M_intern_fifo_rd_en = ((ref_fifo_rd_en == 1'b1) && (first_query_block == 1'b0));

  FIFOWithCount #(
      .DATA_WIDTH(PE_WIDTH),
      .DEPTH_WIDTH(M_FIFO_DEPTH_WIDTH)
  ) M_intern_fifo
  (
      .clk(clk),
      .rst(rst),

      .wr_en(M_intern_fifo_wr_en),
      .rd_en(M_intern_fifo_rd_en),
      .wr_data_in(M_intern_fifo_wr_data_in),

      .rd_data_out(M_intern_fifo_rd_data_out),
      .count(),
      .full(M_intern_fifo_full),
      .empty(M_intern_fifo_empty)
  );

  assign F_intern_fifo_wr_en = init_out[NUM_PE-1];
  assign F_intern_fifo_wr_data_in = F_out[NUM_PE-1];
  assign F_intern_fifo_rd_en = ((ref_fifo_rd_en == 1'b1) && (first_query_block == 1'b0));

  FIFOWithCount #(
      .DATA_WIDTH(PE_WIDTH),
      .DEPTH_WIDTH(F_FIFO_DEPTH_WIDTH)
  ) F_intern_fifo
  (
      .clk(clk),
      .rst(rst),

      .wr_en(F_intern_fifo_wr_en),
      .rd_en(F_intern_fifo_rd_en),
      .wr_data_in(F_intern_fifo_wr_data_in),

      .rd_data_out(F_intern_fifo_rd_data_out),
      .count(),
      .full(F_intern_fifo_full),
      .empty(F_intern_fifo_empty)
  );

  assign bt_logic_start = (next_state == BT_START);
  assign bt_logic_max_score_addr = max_ref_pos_out[NUM_PE-1];
  assign bt_logic_max_score_mod_addr = max_ref_mod_out[NUM_PE-1];
  assign bt_logic_max_score_pe = max_query_pos_out[NUM_PE-1];
  assign bt_logic_max_score_pe_state = max_pe_state_out[NUM_PE-1]; // ????????

  BTLogic #(
      .ADDR_WIDTH(BT_BRAM_ADDR_WIDTH),
      .REF_LEN_WIDTH(REF_BLOCK_SIZE_WIDTH),
      .LOG_NUM_PE(LOG_NUM_PE)
  ) btlogic (
      .clk(clk),
      .rst(rst),
      .start(bt_logic_start),

      .ref_length(bt_ref_length),
      .max_score_addr(bt_logic_max_score_addr),
      .max_score_mod_addr(bt_logic_max_score_mod_addr),
      .max_score_pe(bt_logic_max_score_pe),
      .max_score_pe_state(bt_logic_max_score_pe_state),
      .input_dir(bt_logic_input_dir),
      .input_dir_diag(bt_logic_input_dir_diag),

      .next_addr(bt_logic_next_addr),
      .next_pe(bt_logic_next_pe),
      .next_addr_diag(bt_logic_next_addr_diag),
      .next_pe_diag(bt_logic_next_pe_diag),
      .addr_valid(bt_logic_addr_valid),
      .dir(bt_logic_dir_out),
      .dir_valid(bt_logic_dir_valid),
      .H_offset(bt_logic_H_offset),
      .max_H_offset(max_H_offset),
      .V_offset(bt_logic_V_offset),
      .max_V_offset(max_V_offset),
      .num_tb_steps(bt_logic_num_tb_steps),
      .done(bt_logic_done)
  );

  generate
  for (k = 0; k < NUM_PE; k = k + 1) 
  begin:bram_dir_gen
      assign select_dir_data_in[4*(k+1)-1:4*k] = bt_bram_data_out[k];
      assign bt_bram_addr[k] = (state >= BT_START) ? ((bt_logic_next_pe == k) ? bt_logic_next_addr : bt_logic_next_addr_diag) : pe_dir_addr[k]; // ???
  end
  endgenerate

  mux_1OfN #( 
      .NUM_PORTS_WIDTH(LOG_NUM_PE),
      .DATA_WIDTH(4)
  ) dir_select (
      .clk(clk),
      .select (bt_logic_next_pe), 
      .data_in (select_dir_data_in),
      .data_out (bt_logic_input_dir)
  );


  mux_1OfN #( 
      .NUM_PORTS_WIDTH(LOG_NUM_PE),
      .DATA_WIDTH(4)
  ) dir_select_diag (
      .clk(clk),
      .select (bt_logic_next_pe_diag), 
      .data_in (select_dir_data_in),
      .data_out (bt_logic_input_dir_diag)
  );

  
  always @(*) 
  begin
      next_state = state;
      case (state)
          WAIT:
              if (start) begin
                  next_state = READ_PARAM_FIFO;
              end
          READ_PARAM_FIFO:
              next_state = SET_PARAM;
          SET_PARAM:
              //changing the state from one to last so that the last PE's
              //set_param is changed and the state changes at the same time
              if (set_param[NUM_PE-2] == 1) begin
                  next_state = STREAM_REF_START;
              end
          STREAM_REF_START:
              next_state = STREAM_REF;
          STREAM_REF:
              if (curr_ref_len >= ref_length ) begin
                  next_state = STREAM_REF_STOP;
              end
          STREAM_REF_STOP:
              if (curr_query_len >= query_length) begin
                  next_state = STREAM_REF_DONE;
              end
              else begin
                  next_state = SET_PARAM;
              end
          STREAM_REF_DONE:
              next_state = COMPUTE_MAX_START;
          COMPUTE_MAX_START:
              next_state = COMPUTE_MAX_WAIT;
          COMPUTE_MAX_WAIT:
              if (compute_max_out[NUM_PE-1] == 1'b1) begin
                  if (do_traceback) begin
                      next_state = BT_START;
                  end
                  else begin
                      next_state = DONE;
                  end
              end
          BT_START:
              next_state = BT_WAIT;
          BT_WAIT:
              if (bt_logic_done) begin
                  next_state = DONE;
              end
         DONE:
             next_state = WAIT;
      endcase
  end

endmodule

