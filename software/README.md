## Building 

Use the following steps for compiling on linux machine (tested on Ubuntu 16.04)

* install bond package (https://github.com/Microsoft/bond)
* install tbb package  (latest release can be found at https://github.com/01org/tbb/releases) 
* generate bond headers using gbc (installed with bond library)
```
    $ gbc cpp Darwin.bond
```
* copy kseq.h from https://github.com/lh3/ksw2
* build darwin using the following steps
```
    $ cmake -DCMAKE_BUILD_TYPE=Release .
    $ make
```
* test using the sample data
```
    $ ./Darwin data/sample_ref.fa data/sample_reads.fa  0 > out.maf
```

## Acknowledgment 
Thanks to Roman Snytsar (Microsoft Research) for his contribution to this code.
