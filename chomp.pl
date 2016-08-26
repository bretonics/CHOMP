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
#   FILE:           chomp.pl
#   LICENSE:        MIT
#   USAGE:          Find CRISPR targets and output results for oligo ordering
#   DEPENDENCIES:   - BioPerl modules
#                   - Own 'Modules' repo
#
# ==============================================================================


#-------------------------------------------------------------------------------
# USER VARIABLES
my $DOWNSEQ;
my $UPSEQ;

#-------------------------------------------------------------------------------
# COMMAND LINE
my $SEQ;
my $WINDOWSIZE  = 23;
my $OUTFILE;
my $USAGE       = "\n\n$0 [options]\n
Options:
    -seq                Sequence file to search
    -down               Down sequence to append
    -up                 Up sequence to append
    -window             Window size for CRISPR oligo (default = 23)
    -out                Out file name
    -help               Shows this message
\n";

# OPTIONS
GetOptions(
    'seq=s'             =>\$SEQ,
    'down:s'            =>\$DOWNSEQ,
    'up:s'              =>\$UPSEQ,
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

#-------------------------------------------------------------------------------
# CALLS
my ($CRISPRS, $CRPseqs) = findOligo($sequence, $WINDOWSIZE); #CRISPR HoH and sequences array references
my $CRPfile = writeCRPfasta($CRISPRS, $fileName); #Write CRISPRs FASTA file
my $targets = Search::blast($CRPfile, $SEQ, $WINDOWSIZE); #CRISPR target hits
writeCRPfile($CRISPRS, $targets, $DOWNSEQ, $UPSEQ, $OUTFILE);

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
        die "Did not provide an input file, -seq <infile>", $USAGE;
    }
    unless ($DOWNSEQ){
        say "Did not provide a DOWN stream sequence to append to CRISPR seq, -down <seq>";
    }
    unless ($UPSEQ){
        say "Did not provide an UP stream sequence to append to CRISPR seq, -up <seq>";
    }
    unless ($OUTFILE){
        die "Did not provide an output file, -out <outfile>", $USAGE;
    }
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ();
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes 5 arguments; CRISPR HoH, BLAST target HoA containing matches
# throughout the sequence, the down/up stream sequences to append to crispr
# sequences, and the output file name
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $output = File containing CRISPR target sequences and relative information
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub writeCRPfile {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '(\%CRISPRS, \%targets, $DOWNSEQ, $UPSEQ, $OUTFILE)';
    @_ == 5 or die wrongNumberArguments(), $filledUsage;

    my ($CRISPRS, $targets, $down, $up, $file) = @_;
    my %CRISPRS = %$CRISPRS;
    my %targets = %$targets;
    my $num = keys %CRISPRS; #number of CRISPR sequences

    my $FH = getFH(">", "$OUTDIR/$OUTFILE");
    say $FH "Name\tSequence\tOccurence";

    # Get ordered CRISPR sequences + info to print
    for (my $i = 0; $i < $num; $i++) {
        my $name = "CRISPR_" . $i;
        my $crispr = $CRISPRS{$name}{'oligo'} . $CRISPRS{$name}{'PAM'};
        # Complete oligo sequence:
        # + DOWN flanking target region
        # + CRISPR sequence
        # + UP flanking target region
        my $sequence;
        if (!$down and !$up) {
            $sequence = $crispr;
        } elsif ($down and !$up) {
            $sequence = $down . $crispr;
        } elsif (!$down and $up) {
            $sequence = $crispr . $up;
        } else {
            $sequence = $down . $crispr . $up;
        }


        # CRISPR sequence target matches from BLAST call
        my ($matches) = $targets{$name};
        my $numMatches = @$matches; #number of hashes in array == number of matches for same CRISPR sequence throughout the whole sequence
        say $FH "$name\t$sequence\t$numMatches"; #pring to file
    }
    say "CRISPRs file written to './$OUTDIR/$OUTFILE'";
    return;
}
#-------------------------------------------------------------------------------
# HELPERS

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ($CRISPRS, $OUTDIR, $fileName);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes 3 arguments; HoH reference of CRISPR oligos,
# the output diretory, and the output file name. Writes each CRISPR
# target found in FASTA and returns file location.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = ($outFile);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub writeCRPfasta {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '($CRISPRS, $fileName)';
    @_ == 2 or die wrongNumberArguments(), $filledUsage;

    my ($CRISPRS, $fileName) = @_;
    my $outFile = "$OUTDIR\/$fileName.fasta";
    my $FH = getFH(">", $outFile);
    my $count = 0;

    foreach my $crispr (keys %$CRISPRS) {
        my $oligo = $CRISPRS->{$crispr}->{"oligo"};
        my $PAM = $CRISPRS->{$crispr}->{"PAM"};
        $crispr = $oligo . $PAM ; #join oligo + PAM sequence
        say $FH ">CRISPR_$count\n$crispr";
        $count++;
    } close $FH;

    return $outFile;
}
