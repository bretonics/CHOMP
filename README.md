# CHOMP 🐊


[![GitHub version](https://badge.fury.io/gh/bretonics%2FCHOMP.svg)](http://badge.fury.io/gh/bretonics%2FCHOMP)
[![Github Issues](http://githubbadges.herokuapp.com/bretonics/CHOMP/issues)](https://github.com/bretonics/CHOMP/issues)
[![Pending Pull-Requests](http://githubbadges.herokuapp.com/bretonics/CHOMP/pulls)](https://github.com/bretonics/CHOMP/pulls)
[![GitHub license](https://img.shields.io/badge/License-MIT-red.svg)](https://bretonics.mit-license.org/)
![](https://reposs.herokuapp.com/?path=bretonics/CHOMP&color=lightgrey)



CRISPR Tool
---

CHOMP will search for all **'N'** window sized (default: 23) sequences containing an **NGG** sequence at the tail end, ex.) ATGTAGCTAGCTAGCTAGTA**GGG**.

Will report how many occurrences of this sequence are present in the target sequence (off-target sites), along with base pairs matched for each hit. You can use that to determine which CRISPR target you may want to use.


## Run CHOMP
    perl chomp.pl -seq t/test.fasta -out crisprs


## Arguments
    -seq                Sequence file to search CRISPRs [required]
    -genome             Genome sequence file(s) to BLAST search (search instead of -seq)
    -down               Down sequence to append
    -up                 Up sequence to append
    -window             Window size for CRISPR oligo (default = 23)
    -ss                 Secondary structure prediction
    -out                Out file name [required]
    -outdir             Out directory name
    -html               Write HTML BLAST results to file
    -help               Shows this message


## Output
Writes 2 files under the default directory **./CRISPRS**:

#### **'.FASTA'** file with each CRISPR sequence found


    >CRISPR_0:0
    ATGTAGCTAGCTAGCTAGTAGGG
    >CRISPR_1:23
    AAAAAATTTTCTCTATCTAACGG
    >CRISPR_2:24
    AAAAATTTTCTCTATCTAACGGG
    >CRISPR_3:115
    TGTGATCACGTACTATTATGCGG
    >CRISPR_4:149
    AAAAATCCCATCGATCTAGCAGG
    >CRISPR_5:154
    TCCCATCGATCTAGCAGGCCCGG
    .
    .
    .
    >CRISPR_16:99
    ATAGTACGTGATCACAGTCATGG

Suffix digit after '**:**' denotes nucleotide position in sequence where crispr was found. Ex.) CRISPR_16:**99**
, crispr was found starting at nucleotide position **99** in `-seq` sequence.


#### **'.txt'** report with each CRISPR sequence details

- CHOMP will report how many occurrences of this sequence are present in the target sequence (off-target sites), along with the number of base pair matches (identities) for each. You can use this to determine which CRISPR target is best to use.

| Name | Sequence | Strand | Subject | Start | Occurrences | Identities
| :------------- | :------------- | :------------- | :------------- | :------------- | :------------- | :------------- |
| CRISPR_10 | CTATAGCATGGGCCCCCGATAGG | reverse | test | 230 | 1 | 23
| CRISPR_11 | TATAGCATGGGCCCCCGATAGGG | reverse | test | 229 | 1 | 23
| CRISPR_1 | AAAAAATTTTCTCTATCTAACGG | plus | test | 24 | 3 | 16,10,8
| CRISPR_2 | AAAAATTTTCTCTATCTAACGGG | plus | test | 25 | 3 | 15,10,8
| CRISPR_14 | CTCCGGGCCTGCTAGATCGATGG | reverse | test | 179 | 6 | 15,9,7
| CRISPR_15 | TCCGGGCCTGCTAGATCGATGGG | reverse | test | 178 | 6 | 15,9,7
| CRISPR_5 | TCCCATCGATCTAGCAGGCCCGG | plus | test | 155 | 6 | 15,9,7
| CRISPR_0 | ATGTAGCTAGCTAGCTAGTAGGG | plus | test | 1 | 4 | 14,12,10
| CRISPR_4 | AAAAATCCCATCGATCTAGCAGG | plus | test | 150 | 6 | 11,9,7
| CRISPR_6 | CGAAAAAAAATTTTCCCTATCGG | plus | test | 193 | 1 | 9
| CRISPR_7 | GAAAAAAAATTTTCCCTATCGGG | plus | test | 194 | 1 | 9
| CRISPR_8 | AAAAAAAATTTTCCCTATCGGGG | plus | test | 195 | 1 | 9
| CRISPR_9 | AAAAAAATTTTCCCTATCGGGGG | plus | test | 196 | 1 | 9
| CRISPR_12 | TTCGTCATGCATGCTCGCTCCGG | reverse | test | 196 | 1 | 8
| CRISPR_13 | TCGTCATGCATGCTCGCTCCGGG | reverse | test | 195 | 1 | 8
| CRISPR_16 | ATAGTACGTGATCACAGTCATGG | reverse | test | 131 | 2 | 8
| CRISPR_3 | TGTGATCACGTACTATTATGCGG | plus | test | 116 | 2 | 8

>Table is sorted in decreasing order using top identity for each CRISPR sequence, then sorted by number of occurrences, in current subject.

#### BLAST
Will also write a **blast** results directory under `-outdir` if `-html` option set.
