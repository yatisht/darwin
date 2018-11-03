## Getting started
Run the following commands to perform a reference-guided assembly of long reads (sampled_reads.fa) aligned to a reference (sample_ref.fa) using darwin. The sample reference is chrI of the yeast genome (sacCer3) and sample reads are generated using PBSIM (https://code.google.com/archive/p/pbsim/) for 20X coverage of the reference. Output alignments are in Multiple Alignment Format (MAF).
```
    $ make
    $ ./darwin_ref_guided data/sample_ref.fa data/sample_reads.fa > out.maf
```

To find read overlaps for overlap-layout-consensus based de novo assembly, run following command provided below. Output alignments are in Multiple Alignment Format (MAF). Note that if read A overlaps with read B, only one of A->B or B->A alignment is attempted to be provided.

```
    $ ./read_overlaps data/sample_reads.fa data/sample_reads.fa > out_overlaps.maf
```
