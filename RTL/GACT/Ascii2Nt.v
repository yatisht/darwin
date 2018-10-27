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

// Converting ASCII characters to Nucleotide bases
module Ascii2Nt (
    input [7:0] ascii,
    input complement,
    output reg [3:0] nt);

    localparam A=1, C=2, G=3, T=4, N=0;

    always @(*) begin
        case ({ascii})
            8'h61 : nt = (complement) ? T : A; //8h'61 - ASCII for a 
            8'h41 : nt = (complement) ? T : A; //8h'41 - ASCII for A 
            8'h63 : nt = (complement) ? G : C; //8h'63 - ASCII for c 
            8'h43 : nt = (complement) ? G : C; //8h'43 - ASCII for C 
            8'h67 : nt = (complement) ? C : G; //8h'67 - ASCII for g 
            8'h47 : nt = (complement) ? C : G; //8h'47 - ASCII for G 
            8'h74 : nt = (complement) ? A : T; //8h'74 - ASCII for t 
            8'h54 : nt = (complement) ? A : T; //8h'54 - ASCII for T
            8'h6e : nt = N;
            8'h4e : nt = N;
            default : nt = N;
        endcase
    end
endmodule
