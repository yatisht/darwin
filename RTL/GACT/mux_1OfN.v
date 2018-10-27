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

module mux_1OfN #(
    parameter NUM_PORTS_WIDTH = 2,
    parameter DATA_WIDTH = 32

)(
    input clk,
    input [NUM_PORTS_WIDTH-1:0] select,
    input [(2 ** NUM_PORTS_WIDTH)*DATA_WIDTH-1:0] data_in,
    output reg [DATA_WIDTH-1:0] data_out
);

localparam NUM_PORTS = 2**NUM_PORTS_WIDTH;

always @(posedge clk) begin
    data_out <= data_in[(select+1)*DATA_WIDTH-1-:DATA_WIDTH];
end

endmodule
