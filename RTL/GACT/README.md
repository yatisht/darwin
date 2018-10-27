## RTL Heirarchy
GACTTop.v is the top module. tb_GACTTop.v is the testbench.

## Parameters in testbench
* REF_FILENAME and QUERY_FILENAME consist of the sequences (A, T, C, G, N) in the hex format.
* params consist of the score matrix, gap_open and gap_extend parameters.
* ref_len and query_len are the lengths of the reference and the query sequences respectively.
* score_threshold parameter can be set for filtering out tiles with scores below the threshold.
* GACT can be run in various modes by using the 8 bit align_fields parameter. Its bits represent the following:
  * [7] - don't care
  * [6] - don't care
  * [5] - return direction pointers
  * [4] - reverse reference sequence
  * [3] - complement reference sequence
  * [2] - reverse query sequence
  * [1] - complement query sequence
  * [0] - start at the last cell instead of the maximum score position


## Running test simulation
Running these tests requires the following
* Python
* Icarus Verilog
* GTKWave


This example runs GACT for 10 sets of test data in ./test_data/ . There are 10 reference and query sequences in ./test_data/ref_320.txt and ./test_data/query_320.txt. The 10 refernce and query sequences are sent as in input to GACT in the hex format (./test_data/ref.*.hex, ./test_data/query.*.hex). In order to run GACT simulation, run the following: 
```
    $ ./run_test.sh
```

Results are created in ./results/.
The folder consists of VCD file which consists of the waveform dump. It can be viewed using gtkwave.

```
    $ gtkwave *.vcd
```

The folder also consists of the *.out file which consists of the direction pointer output. It can be further processed to get the aligned sequences, maximum score and maximum score positions. Further processing can be done as follows:

```
    $ python get_alignments.py ./test_data/ref_320.txt ./test_data/query_320.txt 10
```
