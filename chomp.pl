#!/usr/bin/env perl

use strict; use warnings; use diagnostics; use feature qw(say);
use Getopt::Long; use Pod::Usage;

use FindBin; use lib "$FindBin::RealBin/lib";

use Readonly;

use Bio::Seq; use Bio::SeqIO;
use Search;

# Own Modules (https://github.com/bretonics/Modules)
use MyConfig; use MyIO; use Handlers;
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
# COMMAND LINE
my $SEQ;
my @GENOME;
my $DOWNSEQ;
my $UPSEQ;
my $WINDOWSIZE  = 23;
my $OUTFILE;
my $OUTDIR = 'CRISPRS';
my $HTML; #FALSE
my $VERBOSE; #FALSE

my $USAGE       = "\n\n$0 [options]\n
Options:
    -seq                Sequence file to search CRISPRs
    -genome             Genome sequence file to BLAST search (search instead of -seq)
    -down               Down sequence to append
    -up                 Up sequence to append
    -window             Window size for CRISPR oligo (default = 23)
    -out                Out file name
    -outdir             Out directory name
    -html               Print HTML BLAST results
    -help               Shows this message
\n";

# OPTIONS
GetOptions(
    'seq=s'             =>\$SEQ,
    'genome:s{,10}'     =>\@GENOME,
    'down:s'            =>\$DOWNSEQ,
    'up:s'              =>\$UPSEQ,
    'window:i'          =>\$WINDOWSIZE,
    'out=s'             =>\$OUTFILE,
    'outdir:s'          =>\$OUTDIR,
    'html!'             =>\$HTML,
    'verbose!'          =>\$VERBOSE,
    help                =>sub{pod2usage($USAGE);}
)or pod2usage(2);

checks(); #check CL arguments

#-------------------------------------------------------------------------------
# VARIABLES
my $AUTHOR = 'Andres Breton, <dev@andresbreton.com>';

my $REALBIN = $FindBin::RealBin;

my ($seqDetails)    = getSeqDetails($SEQ);
my @SUBJSEQS; # sequence file to use in BLAST search

#-------------------------------------------------------------------------------
# CALLS
mkDir($OUTDIR);
my ($CRISPRS, $CRPseqs) = findOligo($seqDetails, $WINDOWSIZE); # CRISPR HoH and sequences array references
my $CRPfile             = writeCRPfasta($CRISPRS, $OUTFILE); # Write CRISPRs FASTA file
my $targets             = Search::blast($CRPfile, \@SUBJSEQS, $HTML, $OUTFILE, $OUTDIR); # CRISPR target hits
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
        die 'Did not provide an input file, -seq <infile>', $USAGE;
    }
    unless ($DOWNSEQ){
        say 'Did not provide a DOWN stream sequence to append to CRISPR seq, -down <seq>' if($VERBOSE);
    }
    unless ($UPSEQ){
        say 'Did not provide an UP stream sequence to append to CRISPR seq, -up <seq>' if($VERBOSE);
    }
    unless ($OUTFILE){
        die 'Did not provide an output file, -out <outfile>', $USAGE;
    }

    setParameters(); # set parameters to use in calls

    return;
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ();
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes no arguments, uses global variables to decide which set of
# variables to use and set for calls.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $output = sets sequence file to use for BLAST search
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub setParameters {
    # Set up Variables
    if (@GENOME) {
        @SUBJSEQS = @GENOME;
    } elsif ($SEQ) {
        @SUBJSEQS = $SEQ;
    } else {
        die 'Could not determine which file(s) to use as search sequence.'
    }

    return;
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

    my $FH = getFH(">", $outFile);
    say $FH "Name\tSequence\tStrand\tSubject\tStart\tOccurrences\tIdentities";


    my ($sortedCRISPRS, $details) = sortResults(\%targets);
    my %sortedCRISPRS = %$sortedCRISPRS;

    # %details
    # -- Hash key == subject
    # -- Hash key == CRISPR name
    # -- Hash keys == 'identities' and 'occurrences'
    my %details = %$details;

    # Get ordered CRISPR sequences + info to print
    foreach my $crispr (keys %sortedCRISPRS) {
        my @subjects = keys $details{$crispr};
      foreach my $subject ( @subjects) {
        my $sequence = $CRISPRS->{$crispr}->{'sequence'};

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
        my $strand      = $CRISPRS->{$crispr}{'strand'};
        my $occurrence  = $details->{$crispr}{$subject}->{'occurrences'};
        my $identities  = join("," , @{ $details->{$crispr}{$subject}->{'identities'} } ); # get string of identities
        my $sStart      = @{ $targets{$crispr}{$subject} }[0]->{'sstart'}; # get location of BLAST match hit in subject (reference) for CRISPR found
        say $FH "$crispr\t$sequence\t$strand\t$subject\t$sStart\t$occurrence\t$identities"; # print to file
      }
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
# $return = ( \%seqInfo ); Returns sequence details hash
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub getSeqDetails {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '($seqFile)';
    @_ == 1 or die wrongNumberArguments(), $filledUsage;

    my ($seqFile) = @_;

    # Sequence OO
    my $seqInObj    = Bio::SeqIO->new(-file => $seqFile, -alphabet => "dna");
    my $format      = $seqInObj->_guess_format($seqFile); #check format of input file
    my $seqObj      = $seqInObj->next_seq;
    my $sequence    = $seqObj->seq;
    my $reverse     = $seqObj->revcom->seq;

    my %seqInfo = (
        'seqInObj'  => $seqInObj,
        'format'    => $format,
        'seqObj'    => $seqObj,
        'sequence'  => $sequence,
        'reverse'   => $reverse,
    );

    return(\%seqInfo);
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
# ordered array of CRISPR sequence names based on lowest identity
# base pair matches, then occurrences and a details hash with sorted
# identities and number of occurrences per CRISPR sequence.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = ( \@sortedCRISPRS, \%details );
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub sortResults {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '($targetsRef)';
    @_ == 1 or die wrongNumberArguments(), $filledUsage;

    my ($targetsRef) = @_;
    my %targets = %$targetsRef;
    my %details;
    my (%sortedCRISPRS, @sorted);

    # Get number of occurrences per CRISPR target in %targets Hash of Hashes of Array of Hashes
    # -- Hash key == CRISRP name
    # -- Hash key == Subject name
    # -- Array accounts for multiple hits for each CRISPR sequence
    # -- Hash contains BLAST match info

    # CRISPR sequence target matches from BLAST call
    # $targets == (Hash of Hashes of Array of Hashes)
    foreach my $crispr (keys %targets) { # iterate through each CRISPR instance
        my (@subjects) = keys $targets{$crispr}; # @subjects == BLAST subject instances in Hash of Array of Hashes
        foreach my $subject (@subjects) { # iterate through each BLAST subject instance
            my ($matches) = $targets->{$crispr}{$subject};
            my @identities = sortIdentities( $matches );
            # Number of hashes in array == number of matches for same CRISPR sequence throughout the whole sequence
            my $occurrences = @$matches; # number of occcurrences per CRISPR
            # %details
            # -- Hash key == CRISPR name
            # -- Hash key == subject
            # -- Hash keys == 'identities' and 'occurrences'
            $details{$crispr}{$subject} = {
                'identities'    => \@identities,
                'occurrences'   => $occurrences,
            };
        }
    }

    # Return sorted CRISPR names based on lowest identity base pair matches, then occurrences
    foreach my $crispr (keys %details) { # iterate through each CRISPR instance
        my @subjects = keys $details{$crispr};
        # foreach my $subject ( keys $details{$crispr} ) { # iterate through each subject
            @sorted = ( sort { $details{$crispr}{$a}{'identities'}[0] <=> $details{$crispr}{$b}{'identities'}[0] || $details{$crispr}{$a}{'occurrences'} <=> $details{$crispr}{$b}{'occurrences'} } @subjects );
            $sortedCRISPRS{$crispr} = \@sorted;
            # say "$crispr\t" , join ",", $details{$crispr}{$subject}{'identities'}[0];
        # }
        say "$crispr\n", Dumper(\%sortedCRISPRS), "\n\n";
        # foreach (@sorted) {say "$_\n", Dumper($details{$crispr}), "\n\n" }
    }



    return (\%sortedCRISPRS, \%details);
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ($targets{$name});
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes 1 argument: an array reference containing
# all CRISPR matches for a given CRISPR sequence name
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = sorted identities ascending numerically
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub sortIdentities {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '($targets->{$crispr}{$subject})';
    @_ == 1 or die wrongNumberArguments(), $filledUsage;

    my ($matches) = @_;
    my @matches = @$matches;
    my @identities;
    foreach my $hash (@matches) {
        my $nident = $hash->{'nident'}; chomp($nident);
        push @identities, $nident;
    }
    @identities = ( sort {$a <=> $b} @identities ); # sort ascending numerically
    return (@identities); # return sorted identity hits
}
