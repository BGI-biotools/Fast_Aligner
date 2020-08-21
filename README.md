# Fast_Aligner

## Introduction
Large developemnt in Next Generation Sequencing (NGS) technology has enabled DNA sequencing data generating at ever faster speed and at very low cost, which underscores the demand for more efficient tools. Here, we provide **Fast_Aligner**, an efficient and reliable tool, which integrates the pipeline of 'Filteration->Alignment->Sort->Mark Duplicates'.

## Requirements
* htslib 1.9 or greater at [site](https://github.com/samtools/htslib).

## Quick Start Guide
* Build Fast_Aligner: `make`
* Build index for reference genome: `fast_aligner index hg19.fa`
* Align reads to reference: `fast_aligner mem2 hg19.fa input_r1.fq.gz,input_r2.fq.gz`

## Running Fast_Aligner
* **`fast_aligner index`**
This module is designed for building reference index and preparing for the alignment. Besides the path of reference genome fasta file, the `-p` argument can be used to set the output prefix of index files.
* **`fast_aligner mem2`**
The command `mem2` integrates the pipeline of raw reads filteration, reads alignment, sorting alignments and marking PCR duplicates. Besides, this pipeline can be fine-tuned using `-e` argument, which will skip the step of raw reads filteration, and `-x` argument, which will skip the step of markng PCR duplicates.

&ensp;&ensp;&ensp;1. The simplest way to run this module is as follows. Fastq files of mate pairs should be inputted as one parameter and separated using `,`. Please do **not** add any blanks between the two fastq files. Output directory can be set using `-o` argument.
```
fast_aligner mem2 -o out_dir hg19.fa input_r1.fq.gz,input_r2.fq.gz
```

&ensp;&ensp;&ensp;2. Adapter sequences or adapter list files can be added using `-b` argument. Like the fastq files, adapter sequences or adapter list files of mate paris should be inputted as one parameter and separated using `,`, without any blanks between them. Please enter the information of **forward** adapter first, then enter the information of **reverse** adapter.
```
fast_aligner mem2 -o out_dir -b AAGTGACAA,AAGTGAGCCAAGGAGTTG  hg19.fa input_r1.fq.gz,input_r2.fq.gz
```

&ensp;&ensp;&ensp;3. Large sample sequencing data with multiple sequencing libraries or from multiple sequencing lanes can be processed using a `cfg` file.
```
fast_aligner mem2 -o out_dir -i input.cfg hg19.fa
```
&ensp;&ensp;&ensp;Here is an example for `cfg` file.
```
[lib]
  name = 2900978423C
  [lane]
    name = DP800002044BL_L01
    fq1 = DP800002044BL_L01_9_1.fq.gz
    fq2 = DP800002044BL_L01_9_2.fq.gz
    ad1 = AAGTCGGAGGAGACAA
    ad2 = AAGTCGGATGAGCCAAGGAGTTG
  [lane]
    name = DP800002044BL_L02
    fq1 = DP800002044BL_L02_9_1.fq.gz
    fq2 = DP800002044BL_L02_9_2.fq.gz
    ad1 = adatper_1.list.gz
    ad2 = adatper_2.list.gz
    
[lib]
  name = 2900978421C
  [lane]
    name = DP800002047BL_L01
    fq1 = DP800002047BL_L01_11_1.fq.gz
    fq2 = DP800002047BL_L01_11_2.fq.gz
```