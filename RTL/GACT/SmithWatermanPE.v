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

// V <-> H
// F <-> I
// E <-> D
//
module SmithWatermanPE #(
    parameter WIDTH = 10,
    parameter REF_LEN_WIDTH = 10,
    parameter QUERY_LEN_WIDTH = 10,
    parameter LOG_NUM_PE = 2,
    parameter PE_ID = 0
    
)(
    input  clk,                     // System clock
    input  rst,                     // System reset
    
    // Scoring parameters (from substitution matrix) (basically separate indices of the parameter
    // calculated in Nt2Param
    input [WIDTH-1:0] sub_A_in,
    input [WIDTH-1:0] sub_C_in,
    input [WIDTH-1:0] sub_G_in,
    input [WIDTH-1:0] sub_T_in,
    input [WIDTH-1:0] sub_N_in,
    input [WIDTH-1:0] gap_open_in,
    input [WIDTH-1:0] gap_extend_in,
    input set_param,
    
    input  [WIDTH-1:0] V_in,        // Score from previous PE
    input  [WIDTH-1:0] M_in,        // Match score from previous PE
    input  [WIDTH-1:0] F_in,        // Gap penalty of previous PE
    input  [2:0] T_in,              // Reference seq shift in
    input  init_in,                 // Computation active shift in
    input  [WIDTH-1:0] init_V,      // V initialization value
    input  [WIDTH-1:0] init_E,      // E initialization value

    input [REF_LEN_WIDTH-1:0] max_ref_pos_in,
    input [REF_LEN_WIDTH-1:0] max_ref_mod_in,
    input [QUERY_LEN_WIDTH-1:0] max_query_mod_in,
    input [LOG_NUM_PE-1:0] max_query_pos_in,
    input [1:0] max_pe_state_in,
    input compute_max_in,
    input last,
    input last_in,

    output reg [REF_LEN_WIDTH-1:0] max_ref_pos_out,
    output reg [REF_LEN_WIDTH-1:0] max_ref_mod_out,
    output reg [LOG_NUM_PE-1:0] max_query_pos_out,
    output reg [QUERY_LEN_WIDTH-1:0] max_query_mod_out,
    output reg [1:0] max_pe_state_out,

    output reg compute_max_out,
    output reg last_out,
    output wire [WIDTH-1:0] V_out,       // Score of this PE
    output wire [WIDTH-1:0] M_out,       // Match score of this PE
    output wire [WIDTH-1:0] E_out,       // Left gap penalty of this cell
    output wire [WIDTH-1:0] F_out,       // Up Gap penalty of this cell
    output wire [2:0] T_out,             // Reference seq shift out
    output wire init_out,                 // Computation active shift out
    output reg dir_valid,
    output wire [REF_LEN_WIDTH-1:0] dir_addr, 
    output reg signed [3:0] dir 
    );
    
    localparam ZERO=0, MATCH=3, VER=1, HOR=2;

    reg [2:0] T;
    reg signed [WIDTH-1:0] V_diag;
    reg signed [WIDTH-1:0] V;
    reg signed [WIDTH-1:0] M;
    reg signed [WIDTH-1:0] E;
    reg signed [WIDTH-1:0] F;
    reg signed [WIDTH-1:0] max_V;

    reg store_S;
    reg init;
    reg reg_last;

    reg [REF_LEN_WIDTH-1:0] curr_ref_pos;
    reg [REF_LEN_WIDTH-1:0] curr_ref_mod;
    reg [QUERY_LEN_WIDTH-1:0] curr_query_mod;

    reg signed [WIDTH-1:0] sub_A;
    reg signed [WIDTH-1:0] sub_C;
    reg signed [WIDTH-1:0] sub_G;
    reg signed [WIDTH-1:0] sub_T;
    reg signed [WIDTH-1:0] sub_N;
    reg signed [WIDTH-1:0] gap_open;
    reg signed [WIDTH-1:0] gap_extend;
    
    reg signed [WIDTH-1:0] V_gap_open;
    reg signed [WIDTH-1:0] E_gap_extend;
    reg signed [WIDTH-1:0] upV_gap_open;
    reg signed [WIDTH-1:0] upF_gap_extend;
    reg signed [WIDTH-1:0] match_score;
    reg signed [WIDTH-1:0] new_E;
    reg signed [WIDTH-1:0] new_F;
    reg signed [WIDTH-1:0] new_V;
    reg [3:0] new_dir;
    
    reg signed [WIDTH-1:0] match_reward;

    reg [1:0] pe_state;
    reg [1:0] last_pe_state;
    
    assign V_out = V;
    assign M_out = M;
    assign E_out = E;
    assign F_out = F;
    assign T_out = T;
    assign init_out = init;
    assign store_S_out = store_S;
    assign dir_addr = curr_ref_pos-1;

    always @(*) begin
        case ({T_in})//reference sequence base pair for reward calculation (W(r_i,q_j)) to calc match_score
            3'b000 : match_reward = sub_N;
            3'b001 : match_reward = sub_A;
            3'b010 : match_reward = sub_C;
            3'b011 : match_reward = sub_G;
            3'b100 : match_reward = sub_T;
            default : match_reward = 0;
        endcase
    end
    
    //Score calculation
    always @(*) begin
        //for I matrix
        V_gap_open = V + gap_open;
        E_gap_extend = E + gap_extend;
        //for D matrix
        upV_gap_open = V_in + gap_open;
        upF_gap_extend = F_in + gap_extend;
        
        //H(i-1,j-1)+W(r_i,q_j)
        match_score = V_diag + match_reward;          
            
        //I(i,j)
        if ($signed(V_gap_open) > $signed(E_gap_extend)) begin
            new_E = V_gap_open;
            new_dir[3] = 1;
        end
        else begin
            new_E = E_gap_extend;
            new_dir[3] = 0;
        end
        
        //D(i,j)
        if ($signed(upV_gap_open) > $signed(upF_gap_extend)) begin
            new_F = upV_gap_open;
            new_dir[2] = 1;
        end
        else begin
            new_F = upF_gap_extend;
            new_dir[2] = 0;
        end
            
        //calculating max and final state
        if (0 >= $signed(new_E) && 0 >= $signed(new_F) && 0 >= $signed(match_score)) begin
            new_V = 0;
            pe_state = ZERO;
            new_dir[1:0] = 0;
        end
        else if ($signed(new_F) >= $signed(new_E) && $signed(new_F) >= $signed(match_score)) begin
            new_V = new_F;
            pe_state = VER;
            new_dir[1:0] = 1;
        end
        else if ($signed(new_E) >= $signed(match_score)) begin
            new_V = new_E;
            pe_state = HOR;
            new_dir[1:0] = 2;
        end
        else begin
            new_V = match_score;
            pe_state = MATCH;
            new_dir[1:0] = 3;
        end
    end
    //why pe_state and new_dir both? - same value


    always @(posedge clk) begin
        last_pe_state <= pe_state;
        if (rst) begin
            max_ref_pos_out <= 0;
            max_ref_mod_out <= 0;
            max_query_mod_out <= 0;
            max_pe_state_out <= 0;
            max_V <= 0;
        end
        else if ((init == 1'b1) && ((reg_last == 1'b1) || (V > max_V))) begin
            max_ref_pos_out <= curr_ref_pos - 1;
            max_ref_mod_out <= curr_ref_mod - 1;
            max_query_mod_out <= curr_query_mod - 1;
            max_pe_state_out <= last_pe_state;
            max_V <= V;
        end
        else if (compute_max_in) begin
            if (((max_V <= V_in) || (last_in == 1'b1)) && (reg_last == 1'b0)) begin
                max_ref_pos_out <= max_ref_pos_in;
                max_ref_mod_out <= max_ref_mod_in;
                max_query_mod_out <= max_query_mod_in;
                max_pe_state_out <= max_pe_state_in;
            end
        end
    end

    always @(posedge clk) begin
        //reset/initialize variable
        if (rst) begin
            T <= 0;
            V_diag <= 0;
            V <= 0;
            M <= 0;
            E <= (2'b11 << (WIDTH-2));
            F <= 0;
            store_S <= 0;
            init <= 0;
            dir <= 0;
            dir_valid <= 0;
            curr_ref_pos <= 0;
            curr_ref_mod <= 0;
            curr_query_mod <= 0;
            max_query_pos_out <= PE_ID;
            compute_max_out <= 0;
            reg_last <= 0;

        //set the parameters
        end else if (set_param) begin
            sub_A <= sub_A_in;
            sub_C <= sub_C_in;
            sub_G <= sub_G_in;
            sub_T <= sub_T_in;
            sub_N <= sub_N_in;
            gap_open <= gap_open_in;
            gap_extend <= gap_extend_in;
            reg_last <= last;
            init <= 0; 
            dir_valid <= 0;
            curr_ref_mod <= 0;
            curr_query_mod <= curr_query_mod + 1;
        
        //actual calculation
        end else begin
            init <= init_in;
            T <= T_in;
            V_diag <= V_in;
            compute_max_out <= compute_max_in;
            last_out <= (reg_last | last_in);

            if (init_in) begin
                E <= new_E;
                F <= new_F;
                V <= new_V;
                //M <= ($signed(match_score) >= 0) ? match_score : 0;
                dir <= new_dir;
                dir_valid <= 1;
                curr_ref_pos <= curr_ref_pos + 1;
                curr_ref_mod <= curr_ref_mod + 1;
            end else if (compute_max_in) begin
                if ((((max_V > V_in) || (reg_last == 1'b1))) && (last_in == 1'b0)) begin
                    V <= max_V;
                    max_query_pos_out <= PE_ID;
                end
                else begin
                    max_query_pos_out <= max_query_pos_in;
                    V <= V_in;
                end
                dir_valid <= 0;
            end else begin
                V <= init_V;
                E <= init_E;
                dir_valid <= 0;
            end
        end
    end

endmodule

