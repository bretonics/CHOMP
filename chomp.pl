#!/usr/bin/env perl

use strict; use warnings; use diagnostics; use feature qw(say);
use Getopt::Long; use Pod::Usage;

use FindBin; use lib "$FindBin::RealBin/lib";

use Readonly;

use Bio::Seq; use Bio::SeqIO;
use Search;

# Own Modules (https://github.com/bretonics/Modules)
use MyConfig; use MyIO; use Handlers; use Databases;
use Bioinformatics::Eutil;

use Data::Dumper;
# ==============================================================================
#
#   CAPITAN:        Andres Breton, http://andresbreton.com
#   FILE:           chomp.pl :crocodile:
#   LICENSE:        MIT
#   USAGE:          Find CRISPR targets and output results for oligo ordering
#   DEPENDENCIES:   - BioPerl modules
#                   - Own 'Modules' repo
#
# ==============================================================================


#-------------------------------------------------------------------------------
# USER VARIABLES
Readonly my $DW_STREAM = "";
Readonly my $UP_STREAM = "";
#-------------------------------------------------------------------------------
# COMMAND LINE
my $SEQ;
my $WINDOWSIZE  = 23;
my $USAGE       = "\n\n$0 [options]\n
Options:
    -seq                Sequence file to search
    -window             Window size for CRISPR oligo (default = 23)
    -help               Shows this message
\n";

# OPTIONS
GetOptions(
    'seq=s'             =>\$SEQ,
    'window:i'          =>\$WINDOWSIZE,
    help                =>sub{pod2usage($USAGE);}
)or pod2usage(2);
checks(); #check CL arguments
#-------------------------------------------------------------------------------
# VARIABLES
my $AUTHOR = 'Andres Breton, <dev@andresbreton.com>';

my $REALBIN = "$FindBin::RealBin";
my $OUTDIR  = mkDir("CRISPRS");

# Sequence OO
my $seqInObj    = Bio::SeqIO->new(-file => $SEQ, -alphabet => "dna");
my $format      = $seqInObj->_guess_format($SEQ); #check format of input file
my $seqObj      = $seqInObj->next_seq;
my $sequence    = $seqObj->seq;

my ($fileName)  = $SEQ =~ /(\w+)\b\./; #extract file name for output file name
my $seqOutObj   = Bio::SeqIO->new(-file => ">$fileName.fasta", -format => "fasta", -alphabet => "dna");
#-------------------------------------------------------------------------------
# CALLS
my %CRISPRS = findOligo($sequence,$WINDOWSIZE);
# say %CRISPRS;
print Dumper(\%CRISPRS);
#-------------------------------------------------------------------------------
# SUBS

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = checks();
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function checks for arguments passed on the command-line
# using global variables.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = Prompts users and exits if errors
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub checks {
    unless ($SEQ){
        die "\nDid not provide an input file, -file <infile.txt>", $USAGE;
    }
}

#-------------------------------------------------------------------------------
# HELPERS
