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

module Ascii2Param #(
    parameter PE_WIDTH = 10
)
(
    input [7:0] ascii,
    input [12*PE_WIDTH-1:0] in_param,
    output reg [7*PE_WIDTH-1:0] out_param);

    always @(*) begin
        case ({ascii})
            8'h61 : out_param = {in_param[12*PE_WIDTH-1:8*PE_WIDTH], {(PE_WIDTH){1'b0}}, in_param[2*PE_WIDTH-1:0]};
            8'h41 : out_param = {in_param[12*PE_WIDTH-1:8*PE_WIDTH], {(PE_WIDTH){1'b0}}, in_param[2*PE_WIDTH-1:0]};
            8'h63 : out_param = {in_param[11*PE_WIDTH-1:10*PE_WIDTH], in_param[8*PE_WIDTH-1:5*PE_WIDTH], {(PE_WIDTH){1'b0}}, in_param[2*PE_WIDTH-1:0]};
            8'h43 : out_param = {in_param[11*PE_WIDTH-1:10*PE_WIDTH], in_param[8*PE_WIDTH-1:5*PE_WIDTH], {(PE_WIDTH){1'b0}}, in_param[2*PE_WIDTH-1:0]};
            8'h67 : out_param = {in_param[10*PE_WIDTH-1:9*PE_WIDTH], in_param[7*PE_WIDTH-1:6*PE_WIDTH], in_param[5*PE_WIDTH-1:3*PE_WIDTH], {(PE_WIDTH){1'b0}}, in_param[2*PE_WIDTH-1:0]};
            8'h47 : out_param = {in_param[10*PE_WIDTH-1:9*PE_WIDTH], in_param[7*PE_WIDTH-1:6*PE_WIDTH], in_param[5*PE_WIDTH-1:3*PE_WIDTH], {(PE_WIDTH){1'b0}}, in_param[2*PE_WIDTH-1:0]};
            8'h74 : out_param = {in_param[9*PE_WIDTH-1:8*PE_WIDTH], in_param[6*PE_WIDTH-1:5*PE_WIDTH], in_param[4*PE_WIDTH-1:2*PE_WIDTH], {(PE_WIDTH){1'b0}}, in_param[2*PE_WIDTH-1:0]};
            8'h54 : out_param = {in_param[9*PE_WIDTH-1:8*PE_WIDTH], in_param[6*PE_WIDTH-1:5*PE_WIDTH], in_param[4*PE_WIDTH-1:2*PE_WIDTH], {(PE_WIDTH){1'b0}}, in_param[2*PE_WIDTH-1:0]};
            default : out_param = {{(5*PE_WIDTH){1'b0}}, in_param[2*PE_WIDTH-1:0]};
        endcase
    end

endmodule



