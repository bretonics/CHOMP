#!/usr/bin/env perl
return 1 if caller(); #tests

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
#   USAGE:          Find guide RNA (gRNA) sequences for CRISPR targeting
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
my $SS;
my $OUTFILE;
my $VERBOSE;

our $OUTDIR = 'gRNAs';
our $WINDOWSIZE  = 23;


my $USAGE       = "\n\n$0 [options]\n
Options:
    -seq                Sequence file to search gRNAs [required]
    -subjects           Subject sequence file(s) to BLAST search (search instead of -seq)
    -down               Down sequence to append
    -up                 Up sequence to append
    -window             Window size for gRNA sequence (default = 23)
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
my $gRNAs       = findOligo($seqDetails, $WINDOWSIZE); # gRNA HoH
my $fastaFile   = writeFasta($gRNAs, $OUTFILE); # Write gRNAs fasta file
my $targets     = Search::blast($fastaFile, \@SUBJSEQS, $OUTFILE); # gRNA target hits
writeResults($gRNAs, $targets, $DOWNSEQ, $UPSEQ, $WINDOWSIZE, $OUTFILE);

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

    # Exits
    die 'Did not provide an input file, -seq <infile>', $USAGE unless ($SEQ);
    die 'Did not provide an output file, -out <outfile>', $USAGE unless ($OUTFILE) ;

    # Warnings
    say 'Did not provide a DOWN stream sequence to append to gRNA seq, -down <seq>' unless ($DOWNSEQ && $VERBOSE);
    say 'Did not provide an UP stream sequence to append to gRNA seq, -up <seq>' unless ($UPSEQ && $VERBOSE);


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
# This function takes 5 arguments; gRNAs HoH, BLAST target HoA containing matches
# throughout the sequence, the down/up stream sequences to append to gRNA
# sequences, and the output file name
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $output = File containing gRNAs sequences and relative information on results
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub writeResults {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '(\%gRNAs, \%targets, $DOWNSEQ, $UPSEQ, $WINDOWSIZE, $OUTFILE)';
    @_ == 6 or die wrongNumberArguments(), $filledUsage;

    my ($gRNAs, $targets, $down, $up, $window, $file) = @_;

    my %targets = %$targets;
    my $num     = keys %$gRNAs; # number of gRNA sequences
    my $outFile = "$OUTDIR/$file.txt";

    my $FH = getFH(">", $outFile);
    say $FH "Name\tSequence\tStrand\tPalindrome\tSubject\tStart\tOccurrences\tIdentities";

    my ($subjects, $sortedgRNAs, $details) = sortResults(\%targets);
    my @subjects    = @$subjects;
    my %sortedgRNAs = %$sortedgRNAs;
    my %details     = %$details;


    # Get ordered gRNA sequences + info to print
    foreach my $subject (@subjects) {
        foreach my $gRNA ( @{$sortedgRNAs{$subject}} ) {
            my $sequence    = $gRNAs->{$gRNA}->{'sequence'};
            my $palindrome  = $gRNAs->{$gRNA}{'palindrome'};

            # Complete oligo sequence:
            # + DOWN flanking target region
            # + gRNA sequence
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
            my $position    = $gRNAs->{$gRNA}{'start'}; # gRNA sequence position
            my $strand      = $gRNAs->{$gRNA}{'strand'};
            my $occurrence  = $details->{$gRNA}{$subject}->{'occurrences'};
            my $identities  = join("," , @{ $details->{$gRNA}{$subject}->{'unqIdentities'} } ); # get string of identities
            my $sStart      = @{ $targets{$gRNA}{$subject}{'hsps'} }[0]->{'sstart'}; # get location of BLAST match hit in subject (reference) for gRNA found

            say $FH "$gRNA\t$sequence\t$strand\t$palindrome\t$subject\t$sStart\t$occurrence\t$identities"; # print to file
        }
    } close $FH;

    say "gRNAs file written to $outFile";
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
# $input = ($gRNAs, $OUTDIR, $OUTFILE);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes 3 arguments; HoH reference of gRNA sequences,
# the output diretory, and the output file name. Writes each gRNA
# sequence found in FASTA and returns file location.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = ($outFile);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub writeFasta {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '($gRNAs, $OUTFILE)';
    @_ == 2 or die wrongNumberArguments(), $filledUsage;

    my ($gRNAs, $out) = @_;
    my $outFile = "$OUTDIR/$out.fasta";
    my $FH      = getFH(">", $outFile);
    my $num     = keys %$gRNAs; # number of gRNA sequences

    for (my $i = 0; $i < $num; $i++) { # get sequences in numerical order
        my $name = "gRNA_$i";
        my $sequence = $gRNAs->{$name}->{'sequence'};
        ss($name, $sequence) if($SS); # secondary structure prediction if desired
        say $FH ">$name:" . $gRNAs->{$name}->{'start'}; # 'gRNA_0:start_position'
        say $FH $sequence;
    } close $FH;

    return $outFile;
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ($targetsRef);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes 1 arguments: $target hash of array of hashes
# references returned from Search::blast. Returns a numerically
# ordered array of gRNA sequence names based on lowest identity
# base pair matches, then occurrences and a details hash with sorted
# identities and number of occurrences per gRNA sequence.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = ( \@sortedgRNAs, \%details );
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub sortResults {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '($targetsRef)';
    @_ == 1 or die wrongNumberArguments(), $filledUsage;

    my ($targetsRef) = @_;
    # %targets HoHoAoH:
    # -- Hash key == gRNA name
    # -- Hash key == Subject name
    # -- Hash key == 'info'
    # -- Hash key == 'hsps'
    # -- Array accounts for multiple hits for each gRNA sequence
    # -- Hash contains BLAST match info
    my %targets = %$targetsRef;
    my @gRNAs = keys %targets;
    my (@subjects, @noHits, @sorted, %sortedgRNAs, %details);

    # Get number of occurrences from BLAST call per gRNA sequence in %targets Hash of Hashes of Array of Hashes
    foreach my $gRNA (@gRNAs) { # iterate through each gRNA instance
        @subjects = sort keys $targets{$gRNA}; # @subjects == BLAST subject instances in Hash of Arrays of Hash

        foreach my $subject (@subjects) { # iterate through each BLAST subject instance
            # Handle gRNA sequences having 'No hits', skip sorting
            next if ( @{ $targets{$gRNA}{$subject}{'info'} }[0]->{'numhits'} == 0 );

            my @identities;
            my @ids = sortIdentities( $targets{$gRNA}{$subject}{'hsps'} ); # get sorted list of all BLAST hits (and for all subjects) for each gRNA query
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

            # Number of hashes in array == number of matches for same gRNA sequence throughout the whole sequence
            my $occurrences = @identities; # number of occcurrences per gRNA

            # %details
            # -- Hash key == gRNA name
            # -- Hash key == subject
            # -- Hash keys == 'identities', 'unique identities', and 'number of occurrences'
            $details{$gRNA}{$subject} = {
                'identities'    => \@identities,
                'unqIdentities' => \@unqIdentities,
                'occurrences'   => $occurrences,
            };
        }
    }

    # Return sorted gRNA names based on lowest identity base pair matches, then occurrences
    foreach my $subject ( @subjects ) {
        # Get only gRNAs having hits in current subject, remove others
        my @gRNAsHit = grep { @{ $targets{$_}{$subject}{'info'} }[0]->{'numhits'} != 0 } @gRNAs;

        my @sorted = ( sort { $details{$a}{$subject}{'unqIdentities'}[0] <=> $details{$b}{$subject}{'unqIdentities'}[0] || $details{$a}{$subject}{'occurrences'} <=> $details{$b}{$subject}{'occurrences'} } @gRNAsHit );
        $sortedgRNAs{$subject} = \@sorted;
    }

    return (\@subjects, \%sortedgRNAs, \%details);
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ($targets{$name});
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes 1 argument: an array reference containing
# all gRNA matches for a given gRNA sequence name
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = sorted identities ascending numerically
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub sortIdentities {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '($targets{$gRNA}{$subject}{\'hsps\'})';
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
