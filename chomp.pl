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
# COMMAND LINE
my $SEQ;
my $GENOME;
my $DOWNSEQ;
my $UPSEQ;
my $WINDOWSIZE  = 23;
my $OUTFILE;
my $HTML; #FALSE

my $USAGE       = "\n\n$0 [options]\n
Options:
    -seq                Sequence file to search CRISPRs
    -genome             Genome sequence file to BLAST search (search instead of -seq)
    -down               Down sequence to append
    -up                 Up sequence to append
    -window             Window size for CRISPR oligo (default = 23)
    -out                Out file name
    -html               Print HTML BLAST results
    -help               Shows this message
\n";

# OPTIONS
GetOptions(
    'seq=s'             =>\$SEQ,
    'genome:s'          =>\$GENOME,
    'down:s'            =>\$DOWNSEQ,
    'up:s'              =>\$UPSEQ,
    'window:i'          =>\$WINDOWSIZE,
    'out=s'             =>\$OUTFILE,
    'html!'             =>\$HTML,
    help                =>sub{pod2usage($USAGE);}
)or pod2usage(2);

checks(); #check CL arguments

#-------------------------------------------------------------------------------
# VARIABLES
my $AUTHOR = 'Andres Breton, <dev@andresbreton.com>';

my $REALBIN = "$FindBin::RealBin";
my $OUTDIR  = mkDir("crisprs");

my ($seqInObj, $format, $seqObj, $sequence) = setParameters();
my $SEQFILE; #sequence file to use in BLAST search
#-------------------------------------------------------------------------------
# CALLS
my ($CRISPRS, $CRPseqs) = findOligo($sequence, $WINDOWSIZE); #CRISPR HoH and sequences array references
my $CRPfile             = writeCRPfasta($CRISPRS, $OUTFILE); #Write CRISPRs FASTA file
my $targets             = Search::blast($CRPfile, $SEQFILE, $WINDOWSIZE, $HTML); #CRISPR target hits

writeCRPfile($CRISPRS, $targets, $DOWNSEQ, $UPSEQ, $WINDOWSIZE, $OUTFILE);

#-------------------------------------------------------------------------------
# SUBS

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = checks();
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function checks for arguments passed on the command-line using global variables.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = Prompts users and exits if errors
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub checks {
    # Command line arguments passed
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

    setParameters();

    return;
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ();
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes no arguments, uses global variables to decide which set of
# variables to use and set.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $output = sets sequence file to use for BLAST search
# $return = ( $seqInObj, $format, $seqObj, $sequence ); Returns sequence details
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub setParameters {
    # Set up Variables
    if ($GENOME) {
        $SEQFILE = $GENOME;
    } elsif ($SEQ) {
        $SEQFILE = $SEQ;
    } else {
        die "Could not determine wich file to use as search sequence."
    }

    return( getSequence($SEQFILE) );
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
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '(\%CRISPRS, \%targets, $DOWNSEQ, $UPSEQ, $WINDOWSIZE, $OUTFILE)';
    @_ == 6 or die wrongNumberArguments(), $filledUsage;

    my ($CRISPRS, $targets, $down, $up, $window, $file) = @_;
    my %CRISPRS = %$CRISPRS;
    my %targets = %$targets;
    my $num = keys %CRISPRS; #number of CRISPR sequences
    my $outFile = "$OUTDIR/$OUTFILE.txt";

    my $FH = getFH(">", "$outFile");
    say $FH "Name\tSequence\tOccurences\tIdentities(Length:Matches)";

    # Get ordered CRISPR sequences + info to print
    for (my $i = 0; $i < $num; $i++) {
        my $name = "CRISPR_" . $i;
        my $crispr = $CRISPRS{$name}{'gRNA'} . $CRISPRS{$name}{'PAM'};

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

        # Get all percent identities (pident) for each CRISPR match and report
        # how many sequence hits + nucleotide matches for each hit
        my $identities = '';
        foreach my $hash (@$matches) {
            my $pident = $hash->{'pident'}; chomp $pident;
            my $tmp = "$window:$pident";
            $identities = "$identities $tmp,";
        }
        say $FH "$name\t$sequence\t$numMatches\t$identities"; #print to file
    }
    say "CRISPRs file written to './$outFile'";

    return;
}
#-------------------------------------------------------------------------------
# HELPERS

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ($seqFile);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes 1 argument, the sequence file to get details using
# BioPerl modules
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = ( $seqInObj, $format, $seqObj, $sequence ); Returns sequence details
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub getSequence {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '($seqFile)';
    @_ == 1 or die wrongNumberArguments(), $filledUsage;

    my ($seqFile) = @_;

    # Sequence OO
    my $seqInObj    = Bio::SeqIO->new(-file => $seqFile, -alphabet => "dna");
    my $format      = $seqInObj->_guess_format($seqFile); #check format of input file
    my $seqObj      = $seqInObj->next_seq;
    my $sequence    = $seqObj->seq;

    return($seqInObj, $format, $seqObj, $sequence);
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ($CRISPRS, $OUTDIR, $OUTFILE);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes 3 arguments; HoH reference of CRISPR oligos,
# the output diretory, and the output file name. Writes each CRISPR
# target found in FASTA and returns file location.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = ($outFile);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub writeCRPfasta {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '($CRISPRS, $OUTFILE)';
    @_ == 2 or die wrongNumberArguments(), $filledUsage;

    my ($CRISPRS, $OUTFILE) = @_;
    my $outFile = "$OUTDIR\/$OUTFILE.fasta";
    my $FH = getFH(">", $outFile);
    my $count = 0;

    foreach my $crispr (keys %$CRISPRS) {
        my $gRNA = $CRISPRS->{$crispr}->{"gRNA"};
        my $PAM = $CRISPRS->{$crispr}->{"PAM"};
        $crispr = $gRNA . $PAM ; #join gRNA + PAM sequence
        say $FH ">CRISPR_$count\n$crispr";
        $count++;
    } close $FH;

    return $outFile;
}
