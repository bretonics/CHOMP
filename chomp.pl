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

# ==============================================================================
#
#   CAPITAN:        Andres Breton, http://andresbreton.com
#   FILE:           chomp.pl
#   LICENSE:        MIT
#   USAGE:          Find CRISPR targets and output results for oligo ordering
#   DEPENDENCIES:   - BioPerl modules
#                   - Own 'Modules' repo
#
# ==============================================================================


#-------------------------------------------------------------------------------
# USER VARIABLES
Readonly my $DW_STREAM => "DOWNDOWN";
Readonly my $UP_STREAM => "UPUPUP";

#-------------------------------------------------------------------------------
# COMMAND LINE
my $SEQ;
my $WINDOWSIZE  = 23;
my $OUTFILE;
my $USAGE       = "\n\n$0 [options]\n
Options:
    -seq                Sequence file to search
    -window             Window size for CRISPR oligo (default = 23)
    -out                Out file name
    -help               Shows this message
\n";

# OPTIONS
GetOptions(
    'seq=s'             =>\$SEQ,
    'window:i'          =>\$WINDOWSIZE,
    'out=s'             =>\$OUTFILE,
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
my $seqOutObj   = Bio::SeqIO->new(  -file => ">$fileName.fasta",
                                    -format => "fasta",
                                    -alphabet => "dna");

#-------------------------------------------------------------------------------
# CALLS
my $CRISPRS = findOligo($sequence, $WINDOWSIZE); #CRISPR hash of hashes reference
my $CRPfile = writeCRPfasta($CRISPRS, $OUTDIR, $fileName); #CRISPRs FASTA file
my $targets = Search::blast($CRPfile, $SEQ, $WINDOWSIZE); #CRISPR target hits
writeCRPfile($CRISPRS, $targets, $DW_STREAM, $UP_STREAM, $OUTFILE);

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
        die "\nDid not provide an input file, -seq <infile.txt>", $USAGE;
    }
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ();
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes 1 argument;
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $output = File containing CRISPR target information
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub writeCRPfile {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '(\%CRISPRS, \%targets, $DW_STREAM, $UP_STREAM, $OUTFILE)';
    @_ == 5 or die wrongNumberArguments(), $filledUsage;

    my ($CRISPRS, $targets, $down, $up, $file) = @_;
    my @CRPseqs;

    # *************
    # Should move this to findOligo sub in Search module, then return as array + \%CRISPRS
    foreach my $crispr (keys %$CRISPRS)
        push @CRPseqs, $crispr . $CRISPRS->{$crispr}->{"PAM"}; #join oligo + PAM sequence -> push to array
    }

    while ( my ($id, $value) = each(%$targets) ) { #get key-value pair. Value is anonymous array ref of hash(es)
        my $numMatches = @$value; #number of hashes == number of matches for same CRISPR sequence
        # say $down . . $up; next;
        my $FH = getFH(">>", $file);
        while (<$FH>) {
            foreach my $hashRef (@$value) { #each hash is a different match for same CRISPR sequence
                my (%hash) = %$hashRef;

                # say "Crispr is: $id, Hash is: ", $hashRef;
            }
        }
    }

}
#-------------------------------------------------------------------------------
# HELPERS

# ***************** MOVE TO MODULE AS HELPER SUB FOR BLAST
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ($CRISPRS, $OUTDIR, $fileName);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes 3 arguments; HoH reference of CRISPR oligos,
# the output diretory, and the output file name. Writes each CRISPR
# target found in FASTA and returns file location.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = ($outFile);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub writeCRPfasta {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '($CRISPRS, $OUTDIR, $fileName)';
    @_ == 3 or die wrongNumberArguments(), $filledUsage;

    my ($CRISPRS, $OUTDIR, $fileName) = @_;
    my $outFile = "$OUTDIR/$fileName.fasta";
    my $FH = getFH(">", $outFile);
    my $count = 0;

    foreach my $target (keys %$CRISPRS) {
        $target = $target . $CRISPRS->{$target}->{"PAM"}; #join oligo + PAM sequence
        say $FH ">CRISPR_$count\n$target";
        $count++;
    }

    return $outFile;
}
