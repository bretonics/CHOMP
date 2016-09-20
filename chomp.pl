#!/usr/bin/env perl

use strict; use warnings; use diagnostics; use feature qw(say);
use Getopt::Long; use Pod::Usage;

use FindBin; use lib "$FindBin::RealBin/lib";

use Readonly;

use Bio::Seq; use Bio::SeqIO;
use Search;

# Own Modules (https://github.com/bretonics/Modules)
use MyConfig; use MyIO; use Handlers; use Databases;

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
my $OUTDIR  = mkDir("CRISPRS");

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
    my %targets = %$targets;
    my $num = keys %$CRISPRS; #number of CRISPR sequences
    my $outFile = "$OUTDIR/$OUTFILE.txt";

    my $FH = getFH(">", "$outFile");
    say $FH "Name\tSequence\tStrand\tReverse\tOccurrences\tIdentities";


    my ($sortedCRISPRS, $occurrences) = sortOccurrences(\%targets);
    my @sortedCRISPRS = @$sortedCRISPRS;

    # Get ordered CRISPR sequences + info to print
    foreach my $name (@sortedCRISPRS) {
        my $sequence = $CRISPRS->{$name}->{'sequence'};

        # Complete oligo sequence:
        # + DOWN flanking target region
        # + CRISPR sequence
        # + UP flanking target region
        if (!$down and !$up) {
            $sequence = $sequence;
        } elsif ($down and !$up) {
            $sequence = $down . $sequence;
        } elsif (!$down and $up) {
            $sequence = $sequence . $up;
        } else {
            $sequence = $down . $sequence . $up;
        }

        # Get all details to print to file
        my $revSeq = reverse($sequence);
        my $occurrence = $occurrences->{$name};
        my $crispr = $targets{$name}; #Array of Hashes for given CRISPR sequence name
        my $identities = sortIdentities( $crispr );
        my $strand;
        foreach my $hash (@$crispr) {
            my $nident = $hash->{'nident'}; chomp $nident;
            my $tmp = $nident;
            $strand = $hash->{'sstrand'};
        }
        say $FH "$name\t$sequence\t$strand\t$revSeq\t$occurrence\t$identities"; #print to file
    }
    say "CRISPRs file written to $outFile";

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
    my $outFile = "$OUTDIR/$OUTFILE.fasta";
    my $FH = getFH(">", $outFile);
    my $num = keys %$CRISPRS; #number of CRISPR sequences

    for (my $i = 0; $i < $num; $i++) {
        my $crispr = "CRISPR_$i";
        my $sequence = $CRISPRS->{$crispr}->{'sequence'};
        say $FH ">$crispr\n$sequence";
    } close $FH;

    return $outFile;
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ($targetsRef);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes 1 arguments: $target hash of array of hashes
# references returned from Search::blast. Returns a numerically
# ordered array of CRISPR sequence names based on number of
# occurrence and hash with name->occurrence pairing.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = ( \@sortedCRISPRS, \%occurrences );
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub sortOccurrences {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '($targetsRef)';
    @_ == 1 or die wrongNumberArguments(), $filledUsage;

    my ($targetsRef) = @_;
    my %targets = %$targetsRef;
    my %occurrences;
    my @sortedCRISPRS;
    # Get number of occurrences per CRISPR target in %targets Hash of Array of Hashes
    # -- Hash key == CRISRP name
    # -- Array accounts for multiple hits for each CRISPR sequence
    # -- Hash contains BLAST match info
    foreach my $name (keys %targets) {
        # CRISPR sequence target matches from BLAST call
        my ($matches) = $targets{$name}; # $targets == (Hash of Array of Hashes)
        # Number of hashes in array == number of matches for same CRISPR sequence throughout the whole sequence
        $occurrences{$name} = @$matches; # store number of occcurrences per CRISPR
    }

    # Return sorted CRISPR names based on number of occurrences + occurrences
    @sortedCRISPRS = ( sort { $occurrences{$a} <=> $occurrences{$b} } keys %occurrences );
    return (\@sortedCRISPRS, \%occurrences);
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ($targets{$name});
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes 1 argument: an array reference containing
# all CRISPR matches for a given CRISPR sequence name
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = sorted identities descending numerically
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub sortIdentities {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '($targets{$name})';
    @_ == 1 or die wrongNumberArguments(), $filledUsage;

    my ($matches) = @_;
    my @matches = @$matches;
    my @identities;
    foreach my $hash (@matches) {
        my $nident = $hash->{'nident'}; chomp($nident);
        push @identities, $nident;
    }
    @identities = ( sort {$b <=> $a} @identities ); # sort descending numerically
    return join(",", @identities); # return sorted identity hits
}
