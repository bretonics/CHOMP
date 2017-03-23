#!/usr/bin/env perl

use strict; use warnings; use diagnostics; use feature qw(say);

use IO::Detect qw(is_filehandle);
use Test::More;
use Test::Exception; #need this to get dies_ok and lives_ok to work

use FindBin; use lib "$FindBin::RealBin/../lib";

my $SEQ = "usr/test.fasta";
my @GENOME;
my $OUTFILE = "test";;
my $DOWNSEQ = "DWDWDW";
my $UPSEQ = "UPUPUP";
our $OUTDIR = 'tests';
our $WINDOWSIZE  = 23;
my $VERBOSE;
my $SS = "";
our $HTML = "";
my @SUBJSEQS;

mkDir($OUTDIR);

#-------------------------------------------------------------------------------
# TESTS

require_ok 'chomp.pl';

# Modules
BEGIN { use_ok('MyIO') }
BEGIN { use_ok('MyConfig') }
BEGIN { use_ok('Handlers') }
BEGIN { use_ok('Search') }
BEGIN { use_ok('SS') }


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAIN

# # checks
# print "\nTesting checks()...\n";
# dies_ok { checks() } 'dies ok when no arguments passed in checks()';
#
#
# # setParameters
# print "\nTesting setParameters()...\n";
# # dies_ok { setParameters(1) } 'dies ok when 1 arguments is passed in setParameters()';
# dies_ok { setParameters() } 'dies ok when no arguments passed in setParameters()';
# # setParameters();


# getSeqDetails
print "\nTesting getSeqDetails()...\n";
dies_ok { getSeqDetails() } 'dies ok when no argument passed in getSeqDetails()';
dies_ok { getSeqDetails(1) } 'dies ok when 1 incorrect argument is passed in getSeqDetails()';
dies_ok { getSeqDetails(1,2) } 'dies ok when >1 arguments are passed in getSeqDetails()';
dies_ok { getSeqDetails("not.fasta") } 'dies ok when argument is not FASTA file in getSeqDetails()';
lives_ok { getSeqDetails($SEQ) } 'lives ok when file is passed in getSeqDetails()';
my ($seqDetails) = getSeqDetails($SEQ);


# findOligo()
print "\nTesting findOligo()...\n";
dies_ok { findOligo() } 'dies ok when no argument passed in findOligo()';
dies_ok { findOligo(1) } 'dies ok when 1 argument is passed in findOligo()';
dies_ok { findOligo(1, $WINDOWSIZE) } 'dies ok when 2 incorrect arguments passed in findOligo()';
dies_ok { findOligo(1, $WINDOWSIZE) } 'dies ok when 2 arguments passed, but 1 is incorrect in findOligo()';
dies_ok { writeFasta(1,2,3) } 'dies ok when 3 arguments are passed in findOligo()';
lives_ok { findOligo($seqDetails, $WINDOWSIZE) } 'lives ok when right paramaterers passed in findOligo()';

my ($gRNAs, $CRPseqs) = findOligo($seqDetails, $WINDOWSIZE);
is(ref($gRNAs), "HASH", "CRISPRS hash reference returned");
is($gRNAs->{'gRNA_0'}->{'sequence'}, 'ATGTAGCTAGCTAGCTAGTAGGG', 'gRNA_0 sequence is correct');
is($gRNAs->{'gRNA_1'}->{'sequence'}, 'AAAAAATTTTCTCTATCTAACGG', 'gRNA_1 sequence is correct');
is($gRNAs->{'gRNA_2'}->{'sequence'}, 'AAAAATTTTCTCTATCTAACGGG', 'gRNA_2 sequence is correct');
is($gRNAs->{'gRNA_3'}->{'sequence'}, 'TGTGATCACGTACTATTATGCGG', 'gRNA_3 sequence is correct');
is($gRNAs->{'gRNA_4'}->{'sequence'}, 'AAAAATCCCATCGATCTAGCAGG', 'gRNA_4 sequence is correct');
is($gRNAs->{'gRNA_5'}->{'sequence'}, 'TCCCATCGATCTAGCAGGCCCGG', 'gRNA_5 sequence is correct');
is($gRNAs->{'gRNA_6'}->{'sequence'}, 'CGAAAAAAAATTTTCCCTATCGG', 'gRNA_6 sequence is correct');
is($gRNAs->{'gRNA_7'}->{'sequence'}, 'GAAAAAAAATTTTCCCTATCGGG', 'gRNA_7 sequence is correct');
is($gRNAs->{'gRNA_8'}->{'sequence'}, 'AAAAAAAATTTTCCCTATCGGGG', 'gRNA_8 sequence is correct');
is($gRNAs->{'gRNA_9'}->{'sequence'}, 'AAAAAAATTTTCCCTATCGGGGG', 'gRNA_9 sequence is correct');
is($gRNAs->{'gRNA_10'}->{'sequence'}, 'CTATAGCATGGGCCCCCGATAGG', 'gRNA_10 sequence is correct');
is($gRNAs->{'gRNA_11'}->{'sequence'}, 'TATAGCATGGGCCCCCGATAGGG', 'gRNA_11 sequence is correct');
is($gRNAs->{'gRNA_12'}->{'sequence'}, 'TTCGTCATGCATGCTCGCTCCGG', 'gRNA_12 sequence is correct');
is($gRNAs->{'gRNA_13'}->{'sequence'}, 'TCGTCATGCATGCTCGCTCCGGG', 'gRNA_13 sequence is correct');
is($gRNAs->{'gRNA_14'}->{'sequence'}, 'CTCCGGGCCTGCTAGATCGATGG', 'gRNA_14 sequence is correct');
is($gRNAs->{'gRNA_15'}->{'sequence'}, 'TCCGGGCCTGCTAGATCGATGGG', 'gRNA_15 sequence is correct');
is($gRNAs->{'gRNA_16'}->{'sequence'}, 'ATAGTACGTGATCACAGTCATGG', 'gRNA_16 sequence is correct');


# writeFasta
print "\nTesting writeFasta()...\n";
dies_ok { writeFasta() } 'dies ok when no argument passed in writeFasta()';
dies_ok { writeFasta(1) } 'dies ok when 1 argument is passed in writeFasta()';
dies_ok { writeFasta(1,2) } 'dies ok when 2 incorrect arguments passed in writeFasta()';
dies_ok { writeFasta(1, $OUTFILE) } 'dies ok when 2 arguments passed, but 1 is incorrect in writeFasta()';
dies_ok { writeFasta(1,2,3) } 'dies ok when 3 arguments are passed in writeFasta()';
lives_ok { writeFasta($gRNAs, $OUTFILE) } 'lives ok when right parameters are passed in writeFasta()';
my $CRPfile = writeFasta($gRNAs, $OUTFILE);


# Search::blast
print "\nTesting Search::blast()...\n";
dies_ok { Search::blast() } 'dies ok when no argument passed in Search::blast()';
dies_ok { Search::blast(1) } 'dies ok when 1 argument is passed in Search::blast()';
dies_ok { Search::blast(1,2,) } 'dies ok when 2 arguments passed, but at least 1 is incorrect in Search::blast()';
dies_ok { Search::blast(1,2, $OUTFILE) } 'dies ok when 3 arguments passed, but 1 is incorrect in Search::blast()';
dies_ok { Search::blast(1, $CRPfile, $OUTFILE) } 'dies ok when 3 arguments passed, but 2 are incorrect in Search::blast()';
dies_ok { Search::blast(1,2,3) } 'dies ok when 3 incorrect arguments are passed in Search::blast()';
# lives_ok { Search::blast($CRPfile, \@SUBJSEQS, $OUTFILE) } 'lives ok when right parameters are passed in Search::blast()';
# my $targets = Search::blast($CRPfile, \@SUBJSEQS, $OUTFILE);


# _palindrome
print "\nTesting Search::_palindrome()...\n";
dies_ok { Search::_palindrome() } 'dies ok when no argument passed in Search::_palindrome()';
# dies_ok { Search::_palindrome(1) } 'dies ok when argument passed is incorrect in Search::_palindrome()';
lives_ok { Search::_palindrome('ATGTA') } 'lives ok when right parameters are passed in Search::_palindrome()';
is ( Search::_palindrome('ATGTA') , 'Yes', 'YES palindrome found is correct');
is ( Search::_palindrome('ATGTAA') , 'No', 'NO palindrome found is correct');


# writeCRPfile
# print "\nTesting writeCRPfile()...\n";
# writeCRPfile($gRNAs, $targets, $DOWNSEQ, $UPSEQ, $WINDOWSIZE, $OUTFILE);



# Clean Up
`rm -r ./$OUTDIR`;

done_testing();
