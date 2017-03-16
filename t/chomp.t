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

# checks
print "\nTesting checks()...\n";
dies_ok { checks() } 'dies ok when no arguments passed in checks()';


# setParameters
print "\nTesting setParameters()...\n";
# dies_ok { setParameters(1) } 'dies ok when 1 arguments is passed in setParameters()';
dies_ok { setParameters() } 'dies ok when no arguments passed in setParameters()';
# setParameters();


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
dies_ok { writeCRPfasta(1,2,3) } 'dies ok when 3 arguments are passed in findOligo()';
lives_ok { findOligo($seqDetails, $WINDOWSIZE) } 'lives ok when right paramaterers passed in findOligo()';

my ($CRISPRS, $CRPseqs) = findOligo($seqDetails, $WINDOWSIZE);
is(ref($CRISPRS), "HASH", "CRISPRS hash reference returned");
is($CRISPRS->{'CRISPR_0'}->{'sequence'}, 'ATGTAGCTAGCTAGCTAGTAGGG', 'CRISPR_0 sequence is correct');
is($CRISPRS->{'CRISPR_1'}->{'sequence'}, 'AAAAAATTTTCTCTATCTAACGG', 'CRISPR_1 sequence is correct');
is($CRISPRS->{'CRISPR_2'}->{'sequence'}, 'AAAAATTTTCTCTATCTAACGGG', 'CRISPR_2 sequence is correct');
is($CRISPRS->{'CRISPR_3'}->{'sequence'}, 'TGTGATCACGTACTATTATGCGG', 'CRISPR_3 sequence is correct');
is($CRISPRS->{'CRISPR_4'}->{'sequence'}, 'AAAAATCCCATCGATCTAGCAGG', 'CRISPR_4 sequence is correct');
is($CRISPRS->{'CRISPR_5'}->{'sequence'}, 'TCCCATCGATCTAGCAGGCCCGG', 'CRISPR_5 sequence is correct');
is($CRISPRS->{'CRISPR_6'}->{'sequence'}, 'CGAAAAAAAATTTTCCCTATCGG', 'CRISPR_6 sequence is correct');
is($CRISPRS->{'CRISPR_7'}->{'sequence'}, 'GAAAAAAAATTTTCCCTATCGGG', 'CRISPR_7 sequence is correct');
is($CRISPRS->{'CRISPR_8'}->{'sequence'}, 'AAAAAAAATTTTCCCTATCGGGG', 'CRISPR_8 sequence is correct');
is($CRISPRS->{'CRISPR_9'}->{'sequence'}, 'AAAAAAATTTTCCCTATCGGGGG', 'CRISPR_9 sequence is correct');
is($CRISPRS->{'CRISPR_10'}->{'sequence'}, 'CTATAGCATGGGCCCCCGATAGG', 'CRISPR_10 sequence is correct');
is($CRISPRS->{'CRISPR_11'}->{'sequence'}, 'TATAGCATGGGCCCCCGATAGGG', 'CRISPR_11 sequence is correct');
is($CRISPRS->{'CRISPR_12'}->{'sequence'}, 'TTCGTCATGCATGCTCGCTCCGG', 'CRISPR_12 sequence is correct');
is($CRISPRS->{'CRISPR_13'}->{'sequence'}, 'TCGTCATGCATGCTCGCTCCGGG', 'CRISPR_13 sequence is correct');
is($CRISPRS->{'CRISPR_14'}->{'sequence'}, 'CTCCGGGCCTGCTAGATCGATGG', 'CRISPR_14 sequence is correct');
is($CRISPRS->{'CRISPR_15'}->{'sequence'}, 'TCCGGGCCTGCTAGATCGATGGG', 'CRISPR_15 sequence is correct');
is($CRISPRS->{'CRISPR_16'}->{'sequence'}, 'ATAGTACGTGATCACAGTCATGG', 'CRISPR_16 sequence is correct');


# writeCRPfasta
print "\nTesting writeCRPfasta()...\n";
dies_ok { writeCRPfasta() } 'dies ok when no argument passed in writeCRPfasta()';
dies_ok { writeCRPfasta(1) } 'dies ok when 1 argument is passed in writeCRPfasta()';
dies_ok { writeCRPfasta(1,2) } 'dies ok when 2 incorrect arguments passed in writeCRPfasta()';
dies_ok { writeCRPfasta(1, $OUTFILE) } 'dies ok when 2 arguments passed, but 1 is incorrect in writeCRPfasta()';
dies_ok { writeCRPfasta(1,2,3) } 'dies ok when 3 arguments are passed in writeCRPfasta()';
lives_ok { writeCRPfasta($CRISPRS, $OUTFILE) } 'lives ok when right parameters are passed in writeCRPfasta()';
my $CRPfile = writeCRPfasta($CRISPRS, $OUTFILE);


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
# writeCRPfile($CRISPRS, $targets, $DOWNSEQ, $UPSEQ, $WINDOWSIZE, $OUTFILE);



# Clean Up
`rm -r ./$OUTDIR`;

done_testing();
