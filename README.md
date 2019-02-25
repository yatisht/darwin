# darwin
Darwin: A co-processor for long read alignment

## Introduction

This repository provides the software for performing reference-guided or de novo assembly of long reads using algorithms as described in the Darwin manuscript (see reference below). It also provides the RTL code for accelerating the GACT algorithm in hardware (FPGA/ASIC). Refer to the README.md files in software/ and RTL/GACT/ folders to get started.  

## Citing Darwin
* Assembly algorithms and GACT hardware described in:
    
    * Yatish Turakhia, Gill Bejerano, and William J. Dally. "Darwin: A Genomics Co-processor Provides up to 15,000 X Acceleration on Long Read Assembly." Proceedings of the Twenty-Third International Conference on Architectural Support for Programming Languages and Operating Systems (ASPLOS). ACM, 2018. [DOI: 10.1145/3173162.3173193]
    
* Faster seed table construction using AVX instructions in software:
    
    * Roman Snytsar and Yatish Turakhia. "Parallel approach to sliding window sums." arXiv preprint, [arXiv:1811.10074], 2018.
    

## Acknowledgment

* Roman Snytsar ([@mirounga]) for his help restructing the software code to work with the Intel TBB library and for implementing faster seed table construction using AVX instructions.
* George Horrell ([@georgehorrell]) for implementing faster GACT algorithm in software using AVX instructions.
* Sneha D. Goenka ([@gsneha26]) for fixing bugs in GACT RTL.

[arXiv:1811.10074]: https://arxiv.org/abs/1811.10074
[@mirounga]: https://github.com/mirounga
[@georgehorrell]: https://github.com/georgehorrell
[@gsneha26]: https://github.com/gsneha26
