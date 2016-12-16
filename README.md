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

It will report how many occurrences of this sequence are present in the target sequence (off-target sites), along with how many base pair matches for each hit. You can use that to determine which CRISPR target is best to use.


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
Writes 2 files:

1. **.FASTA** file with each CRISPR sequence found
2. **.txt** report with each CRISPR sequence details

Will also write a **blast** directory under `-outdir` if `-html` option set
