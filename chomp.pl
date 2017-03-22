#!/usr/bin/env perl
return 1 if caller();

use strict; use warnings; use diagnostics; use feature qw(say);
use Getopt::Long; use Pod::Usage;
use Readonly;
use Bio::Seq; use Bio::SeqIO;
use FindBin; use lib "$FindBin::RealBin/lib";
use Search; use SS; # lib modules
# Own Modules (https://github.com/bretonics/Modules)
use MyConfig; use MyIO; use Handlers;

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
my @SUBJECTS;
my $DOWNSEQ;
my $UPSEQ;
our $WINDOWSIZE  = 23;
my $SS;
my $OUTFILE;
our $OUTDIR = 'gRNAs';
my $VERBOSE;

my $USAGE       = "\n\n$0 [options]\n
Options:
    -seq                Sequence file to search CRISPRs [required]
    -subjects           Subject sequence file(s) to BLAST search (search instead of -seq)
    -down               Down sequence to append
    -up                 Up sequence to append
    -window             Window size for CRISPR oligo (default = 23)
    -ss                 Secondary structure prediction
    -out                Out file name [required]
    -outdir             Out directory name
    -help               Shows this message
\n";

# OPTIONS
GetOptions(
    'seq=s'             =>\$SEQ,
    'subjects:s{,10}'   =>\@SUBJECTS,
    'down:s'            =>\$DOWNSEQ,
    'up:s'              =>\$UPSEQ,
    'window:i'          =>\$WINDOWSIZE,
    'ss!'               =>\$SS,
    'out=s'             =>\$OUTFILE,
    'outdir:s'          =>\$OUTDIR,
    'verbose!'          =>\$VERBOSE,
    help                =>sub{pod2usage($USAGE);}
)or pod2usage(2);

checks(); # check CL arguments

#-------------------------------------------------------------------------------
# VARIABLES
my $AUTHOR = 'Andres Breton, <dev@andresbreton.com>';

my $REALBIN = $FindBin::RealBin;

my ($seqDetails) = getSeqDetails($SEQ);
my @SUBJSEQS; # sequence file to use in BLAST search

#-------------------------------------------------------------------------------
# CALLS
mkDir($OUTDIR);
my $CRISPRS             = findOligo($seqDetails, $WINDOWSIZE); # CRISPR HoH
my $CRPfile             = writeCRPfasta($CRISPRS, $OUTFILE); # Write CRISPRs FASTA file
my $targets             = Search::blast($CRPfile, \@SUBJSEQS, $OUTFILE); # CRISPR target hits
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
    unless ($SEQ) {
        die 'Did not provide an input file, -seq <infile>', $USAGE;
    }
    unless ($DOWNSEQ){
        say 'Did not provide a DOWN stream sequence to append to CRISPR seq, -down <seq>' if($VERBOSE);
    }
    unless ($UPSEQ){
        say 'Did not provide an UP stream sequence to append to CRISPR seq, -up <seq>' if($VERBOSE);
    }
    unless ($OUTFILE) {
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
    if (@SUBJECTS) {
        @SUBJSEQS = @SUBJECTS;
    } elsif ($SEQ) {
        @SUBJSEQS = $SEQ;
    } else {
        die 'Could not determine which file(s) to use as search sequence.'
    }

    mkDir($OUTDIR); mkDir("$OUTDIR/ss") if($SS);

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
    my $num     = keys %$CRISPRS; # number of CRISPR sequences
    my $outFile = "$OUTDIR/$OUTFILE.txt";

    my $FH = getFH(">", $outFile);
    say $FH "Name\tSequence\tStrand\tPalindrome\tSubject\tStart\tOccurrences\tIdentities";

    my ($subjects, $sortedCRISPRS, $details) = sortResults(\%targets);
    my @subjects        = @$subjects;
    my %sortedCRISPRS   = %$sortedCRISPRS;
    my %details         = %$details;


    # Get ordered CRISPR sequences + info to print
    foreach my $subject (@subjects) {
        foreach my $crispr ( @{$sortedCRISPRS{$subject}} ) {
            my $sequence = $CRISPRS->{$crispr}->{'sequence'};
            my $palindrome = $CRISPRS->{$crispr}{'palindrome'};

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
            my $position    = $CRISPRS->{$crispr}{'start'}; #CRISPR sequence position
            my $strand      = $CRISPRS->{$crispr}{'strand'};
            my $occurrence  = $details->{$crispr}{$subject}->{'occurrences'};
            my $identities  = join("," , @{ $details->{$crispr}{$subject}->{'unqIdentities'} } ); # get string of identities
            my $sStart      = @{ $targets{$crispr}{$subject}{'hsps'} }[0]->{'sstart'}; # get location of BLAST match hit in subject (reference) for CRISPR found

            say $FH "$crispr\t$sequence\t$strand\t$palindrome\t$subject\t$sStart\t$occurrence\t$identities"; # print to file
        }
    } close $FH;

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
    my $sequence    = uc $seqObj->seq;
    my $reverse     = uc $seqObj->revcom->seq;

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
    my $num = keys %$CRISPRS; # number of CRISPR sequences

    for (my $i = 0; $i < $num; $i++) { # get sequences in numerical order
        my $crispr = "gRNA_$i";
        my $sequence = $CRISPRS->{$crispr}->{'sequence'};
        ss($crispr, $sequence) if($SS); # secondary structure prediction if desired
        say $FH ">$crispr:" . $CRISPRS->{$crispr}->{'start'}; # 'CRISPR_0:start_position'
        say $FH $sequence;
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
    # %targets HoHoAoH:
    # -- Hash key == CRISRP name
    # -- Hash key == Subject name
    # -- Hash key == 'info'
    # -- Hash key == 'hsps'
    # -- Array accounts for multiple hits for each CRISPR sequence
    # -- Hash contains BLAST match info
    my %targets = %$targetsRef;
    my @crisprs = keys %targets;
    my (@subjects, @noHits, @sorted, %sortedCRISPRS, %details);

    # Get number of occurrences from BLAST call per CRISPR target in %targets Hash of Hashes of Array of Hashes
    foreach my $crispr (@crisprs) { # iterate through each CRISPR instance
        @subjects = sort keys $targets{$crispr}; # @subjects == BLAST subject instances in Hash of Arrays of Hash

        foreach my $subject (@subjects) { # iterate through each BLAST subject instance
            # Handle CRISPR sequences having 'No hits', skip sorting
            next if ( @{ $targets{$crispr}{$subject}{'info'} }[0]->{'numhits'} == 0 );

            my @identities;
            my @ids = sortIdentities( $targets->{$crispr}{$subject}{'hsps'} ); # get sorted list of all BLAST hits (and for all subjects) for each CRISPR query
            push @identities, @ids; # push identities list for each subject

            # Remove duplicates
            my @unqIdentities;
            foreach my $value (@identities) {
                # Push identitities value unless already recorded
                next if ( grep { $_ == $value} @unqIdentities );
                push @unqIdentities, $value;
            }

            # Sort identities from all HSPs
            @identities     = ( sort {$b <=> $a} @identities );
            @unqIdentities  = ( sort {$b <=> $a} @unqIdentities );

            # Number of hashes in array == number of matches for same CRISPR sequence throughout the whole sequence
            my $occurrences = @identities; # number of occcurrences per CRISPR

            # %details
            # -- Hash key == CRISPR name
            # -- Hash key == subject
            # -- Hash keys == 'identities', 'unique identities', and 'number of occurrences'
            $details{$crispr}{$subject} = {
                'identities'    => \@identities,
                'unqIdentities' => \@unqIdentities,
                'occurrences'   => $occurrences,
            };
        }
    }

    # Return sorted CRISPR names based on lowest identity base pair matches, then occurrences
    foreach my $subject ( @subjects ) {
        # Get only CRISPRs having hits in current subject, remove others
        my @crisprsHit = grep { @{ $targets{$_}{$subject}{'info'} }[0]->{'numhits'} != 0 } @crisprs;

        my @sorted = ( sort { $details{$a}{$subject}{'unqIdentities'}[0] <=> $details{$b}{$subject}{'unqIdentities'}[0] || $details{$a}{$subject}{'occurrences'} <=> $details{$b}{$subject}{'occurrences'} } @crisprsHit );
        $sortedCRISPRS{$subject} = \@sorted;
    }

    return (\@subjects, \%sortedCRISPRS, \%details);
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
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '($targets->{$crispr}{$subject}{\'hsps\'})';
    @_ == 1 or die wrongNumberArguments(), $filledUsage;

    my ($hsps) = @_;
    my @hsps = @$hsps;
    my $numMatches = @hsps;
    my @identities;

    foreach my $hash (@hsps) {
        my $nident = $hash->{'nident'}; chomp($nident);
        return ($nident) if ($numMatches == 1); # return single match (when $nident == $WINDOWSIZE is only match)
        push @identities, $nident unless ($nident == $WINDOWSIZE);
    }

    @identities = ( sort {$b <=> $a} @identities ); # sort descending numerically
    return (@identities); # return sorted identity hits
}
