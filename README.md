#CHOMP 🐊


[![GitHub version](https://badge.fury.io/gh/bretonics%2FCHOMP.svg)](http://badge.fury.io/gh/bretonics%2FCHOMP)
[![Github Issues](http://githubbadges.herokuapp.com/bretonics/CHOMP/issues)](https://github.com/bretonics/CHOMP/issues)
[![Pending Pull-Requests](http://githubbadges.herokuapp.com/bretonics/CHOMP/pulls)](https://github.com/bretonics/CHOMP/pulls)
[![GitHub license](https://img.shields.io/badge/License-MIT-red.svg)](https://bretonics.mit-license.org/)
![](https://reposs.herokuapp.com/?path=bretonics/CHOMP&color=lightgrey)



CRISPR Tool
---
Search for CRISPR sequences given a genome or gene sequence.

CHOMP will search for all **'N'** window sized (default: 23) sequences containing an **NGG** sequence at the tail end, ex.) ATGTAGCTAGCTAGCTAGTA**GGG**.

Will report how many occurrences of this sequence are present in the target sequence (off-target sites), along with how many base pair matches for each hit. You can use that to determine which CRISPR target you may want to use.


##Arguments
	-seq                Sequence file to search CRISPRs
    -genome             Genome sequence file to BLAST search (search instead of -seq)
    -down               Down sequence to append
    -up                 Up sequence to append
    -window             Window size for CRISPR oligo (default = 23)
    -out                Out file name
    -html               Print HTML BLAST results
    -help               Shows this message
    

##Output
Writes 2 files:

1. FASTA file with each CRISPR sequence found
2. Report with each CRISPR sequence details