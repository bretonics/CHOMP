# CHOMP 🐊

[![GitHub version](https://badge.fury.io/gh/bretonics%2FCHOMP.svg)](https://badge.fury.io/gh/bretonics%2FCHOMP)
[![GitHub license](https://img.shields.io/badge/License-MIT-red.svg)](https://bretonics.mit-license.org/)
[![codecov](https://codecov.io/gh/bretonics/CHOMP/branch/unit-tests/graph/badge.svg)](https://codecov.io/gh/bretonics/CHOMP/branch/unit-tests)
[![Github Issues](https://githubbadges.herokuapp.com/bretonics/CHOMP/issues)](https://github.com/bretonics/CHOMP/issues)
[![Pending Pull-Requests](https://githubbadges.herokuapp.com/bretonics/CHOMP/pulls)](https://github.com/bretonics/CHOMP/pulls)
![](https://reposs.herokuapp.com/?path=bretonics/CHOMP&color=lightgrey)

CRISPR Tool
---

CHOMP will search for all **'n'** window sized guide RNAs (gRNAs) sequences containing an **NGG**  at the tail end. Default window size is 23.

> ex.) ATGTAGCTAGCTAGCTAGTA**GGG**.

It will report how many occurrences of this sequence are present in the target sequence (off-target sites), along with number of matching bases for each subject/query hit. You can use that to determine which gRNAs you may want to use as CRISPR target(s).

## Run CHOMP

    perl chomp.pl -seq usr/test.fasta -out gRNAs

## Arguments

    -seq                Sequence file to search gRNAs [required]
    -genome             Genome sequence file(s) to BLAST search (search instead of -seq)
    -down               Down sequence to append
    -up                 Up sequence to append
    -window             Window size for gRNA oligo (default = 23)
    -ss                 Secondary structure prediction
    -out                Out file name [required]
    -outdir             Out directory name
    -help               Shows this message

## Output

Writes 2 files under the default directory **gRNAs**:

#### **'.fasta'**

Fasta file of each gRNA sequence found.

    >gRNA_0:0
    ATGTAGCTAGCTAGCTAGTAGGG
    >gRNA_1:23
    AAAAAATTTTCTCTATCTAACGG
    >gRNA_2:24
    AAAAATTTTCTCTATCTAACGGG
    >gRNA_3:115
    TGTGATCACGTACTATTATGCGG
    >gRNA_4:149
    AAAAATCCCATCGATCTAGCAGG
    >gRNA_5:154
    TCCCATCGATCTAGCAGGCCCGG
    .
    .
    .
    >gRNA_16:99
    ATAGTACGTGATCACAGTCATGG

Suffix digit after '**:**' denotes nucleotide position in sequence where gRNA was found. Ex.) gRNA_16:**99**
, gRNA was found starting at nucleotide position **99** in `-seq` sequence.

#### **'.txt'**

Report with each gRNA sequence's details.

CHOMP will report how many occurrences of this sequence are present in the target sequence (off-target sites), along with the number of base pair matches (identities) for each. You can use this to determine which gRNA sequence is best to use for target.

| Name | Sequence | Strand | Palindrome | Subject | Start | Occurrences | Identities |
| :------------- | :------------- | :------------- | :------------- | :------------- | :------------- | :------------- | :------------- |
| gRNA_13 | TCGTCATGCATGCTCGCTCCGGG | reverse | No | test | 173 | 1 | 8
| gRNA_12 | TTCGTCATGCATGCTCGCTCCGG | reverse | No | test | 174 | 1 | 8
| gRNA_3 | TGTGATCACGTACTATTATGCGG | plus | No | test | 116 | 2 | 8
| gRNA_16 | ATAGTACGTGATCACAGTCATGG | reverse | No | test | 109 | 2 | 8
| gRNA_8 | AAAAAAAATTTTCCCTATCGGGG | plus | No | test | 195 | 1 | 9
| gRNA_7 | GAAAAAAAATTTTCCCTATCGGG | plus | No | test | 194 | 1 | 9
| gRNA_6 | CGAAAAAAAATTTTCCCTATCGG | plus | No | test | 193 | 1 | 9
| gRNA_9 | AAAAAAATTTTCCCTATCGGGGG | plus | No | test | 196 | 1 | 9
| gRNA_4 | AAAAATCCCATCGATCTAGCAGG | plus | No | test | 150 | 6 | 11,9,7
| gRNA_0 | ATGTAGCTAGCTAGCTAGTAGGG | plus | No | test | 1 | 4 | 14,12,10
| gRNA_2 | AAAAATTTTCTCTATCTAACGGG | plus | No | test | 25 | 3 | 15,10,8
| gRNA_15 | TCCGGGCCTGCTAGATCGATGGG | reverse | No | test | 156 | 6 | 15,9,7
| gRNA_14 | CTCCGGGCCTGCTAGATCGATGG | reverse | No | test | 157 | 6 | 15,9,7
| gRNA_5 | TCCCATCGATCTAGCAGGCCCGG | plus | No | test | 155 | 6 | 15,9,7
| gRNA_1 | AAAAAATTTTCTCTATCTAACGG | plus | No | test | 24 | 3 | 16,10,8
| gRNA_11 | TATAGCATGGGCCCCCGATAGGG | reverse | No | test | 207 | 1 | 23
| gRNA_10 | CTATAGCATGGGCCCCCGATAGG | reverse | No | test | 208 | 1 | 23

>Table is sorted in increasing order using the top identity for each gRNA sequence, and then sorted by number of occurrences, in current subject.
