#!/usr/bin/env perl

use strict; use warnings; use diagnostics; use feature qw(say);
use Getopt::Long; use Pod::Usage;

use FindBin; use lib "$FindBin::RealBin/lib";

use Readonly;

use Bio::Seq; use Bio::SeqIO;

# Own Modules (https://github.com/bretonics/Modules)
use Bioinformatics::MyConfig;
use Bioinformatics::MyIO;
use Bioinformatics::Eutil;
use Databases; use MyConfig;


# ==============================================================================
#
#   CAPITAN: Andres Breton http://andresbreton.com
#   FILE: chomp.pl 🐊
#   LICENSE: MIT
#   USAGE: Find CRISPR targets and output results for oligo ordering
#   DEPENDENCIES:   - BioPerl modules
#                   - Own Modules git repo
#
# ==============================================================================


#-------------------------------------------------------------------------------
# USER VARIABLES
Readonly my $DW_STREAM = "";
Readonly my $UP_STREAM = "";
#-------------------------------------------------------------------------------
# COMMAND LINE
my $SEQ;
my @IDS;
my $MONGODB = "CRISPR";
my $COLLECTION;
my ($INSERT, @UPDATE, @READ, @REMOVE);
my $USAGE= "\n\n$0 [options]\n
Options:
    -seq                Sequence file to search
    -help               Shows this message
\n";

# OPTIONS
GetOptions(
    'seq=s'             =>\$SEQ,
    help                =>sub{pod2usage($USAGE);}
)or pod2usage(2);
checks(); #check CL arguments
#-------------------------------------------------------------------------------
# VARIABLES
my $REALBIN = "$FindBin::RealBin";
my $OUTDIR = mkDir("CRISPRS");

# Color Output
my $GRNTXT = "\e[1;32m"; #bold green
my $REDTXT = "\e[1;31m"; #bold red
my $NC = "\e[0m"; #color reset

# Sequence OO
# my $seqInObject = Bio::SeqIO->new(-file => $SEQ, -format => "genbank", -alphabet => "dna");
# my $format = $seqInObject->_guess_format($SEQ); #check format of input file
#
# my $sequence = $seqInObject->next_seq;
# my $actual = $sequence->seq;
#
# my ($fileName) = $SEQ =~ /(\w+)\b\./; #extract file name
# my $seqOutObject = Bio::SeqIO->new(-file => ">$fileName.fasta", -format => "fasta", -alphabet => "dna");
#-------------------------------------------------------------------------------
# CALLS

#-------------------------------------------------------------------------------
# SUBS
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = checks();
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function checks for arguments passed on the command-line
# using global variables.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = Prompts users and exits if errors
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub checks {
    unless ($SEQ){
        die "Did not provide an input file, -file <infile.txt>", $USAGE;
    }
    if ($format ne "genbank") {
        say "Incorrect file format. You provided $format format.";
        say "Please provide a GenBank file format";
    }

}

#-------------------------------------------------------------------------------
# HELPERS
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ("DirName");
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes one argument, a string to write directory
# if non-existent.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = ();
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub _mkDir {
    my ($outDir) = @_;
    `mkdir $outDir` unless(-e $outDir);
}
