# CHOMP 🐊


[![GitHub version](https://badge.fury.io/gh/bretonics%2FCHOMP.svg)](http://badge.fury.io/gh/bretonics%2FCHOMP)
[![Github Issues](http://githubbadges.herokuapp.com/bretonics/CHOMP/issues)](https://github.com/bretonics/CHOMP/issues)
[![Pending Pull-Requests](http://githubbadges.herokuapp.com/bretonics/CHOMP/pulls)](https://github.com/bretonics/CHOMP/pulls)
[![GitHub license](https://img.shields.io/badge/License-MIT-red.svg)](https://bretonics.mit-license.org/)
![](https://reposs.herokuapp.com/?path=bretonics/CHOMP&color=lightgrey)



CRISPR Tool
---
>Search for all CRISPR sequences provided a genome or gene sequence.

CHOMP will search for all **'N'** window sized (default: 23) sequences containing an **NGG** sequence at the tail end, ex.) ATGTAGCTAGCTAGCTAGTA**GGG**.

Will report how many occurrences of this sequence are present in the target sequence (off-target sites), along with how many base pair matches for each hit. You can use that to determine which CRISPR target you may want to use.


## Arguments
    -seq                Sequence file to search CRISPRs
    -genome             Genome sequence file to BLAST search (search instead of -seq)
    -down               Down sequence to append
    -up                 Up sequence to append
    -window             Window size for CRISPR oligo (default = 23)
    -out                Out file name
    -html               Print HTML BLAST results
    -help               Shows this message
    
## Run CHOMP
    perl chomp.pl -seq t/test.fasta -out crisprs
    
    
## Output
Writes 2 files under the default directory **./CRISPRS**:

#### FASTA file with each CRISPR sequence found


    >CRISPR_0
    ATGTAGCTAGCTAGCTAGTAGGG
    >CRISPR_1
    AAAAAATTTTCTCTATCTAACGG
    >CRISPR_2
    AAAAATTTTCTCTATCTAACGGG
    >CRISPR_3
    TGTGATCACGTACTATTATGCGG
    >CRISPR_4
    AAAAATCCCATCGATCTAGCAGG
    >CRISPR_5
    TCCCATCGATCTAGCAGGCCCGG
    >CRISPR_6
    CGAAAAAAAATTTTCCCTATCGG
    >CRISPR_7
    GAAAAAAAATTTTCCCTATCGGG
    >CRISPR_8
    AAAAAAAATTTTCCCTATCGGGG
    >CRISPR_9
    AAAAAAATTTTCCCTATCGGGGG

#### Report with each CRISPR sequence details

| Name            | Sequence        | SubjStart       | Occurrences     | Identities
| :------------- | :------------- | :------------- | :------------- | :------------- |
| CRISPR_8 | AAAAAAAATTTTCCCTATCGGGG |  195 | 2 | 23,9
| CRISPR_9 | AAAAAAATTTTCCCTATCGGGGG |  196 | 2 | 23,9
| CRISPR_6 | CGAAAAAAAATTTTCCCTATCGG |  193 | 2 | 23,9
| CRISPR_7 | GAAAAAAAATTTTCCCTATCGGG |  194 | 2 | 23,9
| CRISPR_3 | TGTGATCACGTACTATTATGCGG |  116 | 3 | 23,8,8
| CRISPR_2 | AAAAATTTTCTCTATCTAACGGG |  25  | 4 | 23,15,8,8
| CRISPR_1 | AAAAAATTTTCTCTATCTAACGG |  24  | 4 | 23,16,8,8
| CRISPR_0 | ATGTAGCTAGCTAGCTAGTAGGG |  1   | 5 | 23,14,12,10,10
| CRISPR_5 | TCCCATCGATCTAGCAGGCCCGG |  155 | 7 | 23,15,9,7,7,7,7
| CRISPR_4 | AAAAATCCCATCGATCTAGCAGG |  150 | 8 | 23,9,7,7,7,7,7,7
