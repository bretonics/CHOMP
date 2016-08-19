package Search;

use Exporter qw(import);
our @ISA = qw(Exporter);
our @EXPORT = qw(findOligo); #functions exported by default
our @EXPORT_OK = qw(blast); #functions for explicit export

use strict; use warnings; use diagnostics; use feature qw(say);
use Carp;

use MyConfig; use MyIO;

# ==============================================================================
#
#   CAPITAN:        Andres Breton, http://andresbreton.com
#   FILE:           search.pm
#   LICENSE:
#   USAGE:
#   DEPENDENCIES:   - NCBI's BLAST+ CL utility
#
# ==============================================================================

=head1 NAME

Search - package searching CRISPR sequences and offsite targets

=head1 SYNOPSIS

Creation:
    use Search;

=head1 DESCRIPTION


=head1 EXPORTS

=head2 Default Behaviors

Exports findOligo subroutine by default

use Search;

=head2 Optional Behaviors

Search::blast;

=head1 FUNCTIONS

=cut

#-------------------------------------------------------------------------------
# MAIN

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
=head2 findOligo

    Arg [1]     : Sequence to be searched

    Arg [2]     : Window size of CRISPR target

    Example     : findOligo($sequence, $windowSize)

    Description : Find CRISPR targets

    Returntype  : Hash of hashes reference

    Status      : Development

=cut
sub findOligo {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '($sequence, $windowSize)';
    @_ == 2 or confess wrongNumberArguments(), $filledUsage;

    my ($sequence, $windowSize) = @_;
    my $seqLen = length($sequence);

    my (%CRISPRS, @CRPseqs, $oligo, $PAM, $content, $contentG, $contentC, $GC);
    my $instance = 0; #track CRISPR count

    for (my $i = 0; $i < $seqLen; $i++) {
        my $window = substr $sequence, $i, $windowSize;

        # Return CRISPR sequences and information once done
        if ( length($window) < $windowSize ) { #don't go out of bounds when at end of sequence, return CRISPR sequences found
            foreach my $name (keys %CRISPRS) {
                my $crispr = $CRISPRS{$name}{"oligo"} . $CRISPRS{$name}{"PAM"}; #join oligo + PAM sequence
                push @CRPseqs, $crispr #push to array
            }
            # Return references of HoH containing all CRISPR instances found and respective information for each and array with just the sequences joined (kmer oligo + PAM)
            return(\%CRISPRS, \@CRPseqs);
        };

        my $kmer = ($windowSize - 3); #kmer is the string of base pairs before NGG

        if ($window =~ /(.+)(.GG)$/) {
            ($oligo, $PAM) = ($1, $2); #get first 'kmer' number of nucleotides in oligo + PAM (NGG)
            my $name = "CRISPR_$instance"; $instance++;

            # GC Content
            $contentG = $window =~ tr/G//;
            $contentC = $window =~ tr/C//;
            $GC = ($contentG + $contentC)/$windowSize;

            # Store CRISPR oligomers and info in Hash of Hashes
            $content = { #anonymous hash of relevant oligo content
                'oligo' => $oligo,
                'PAM'   => $PAM,
                'G'     => $contentG,
                'C'     => $contentC,
                'GC'    => $GC,
            };
            # Hash key == CRISPR sequence name
            # Hash value == HoH with CRISPR content info
            $CRISPRS{$name} = $content;
        }
    }
}
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
=head2 blast

    Arg [1]     : CRISPR hash reference returned in findOligo sub

    Arg [2]     : Sequence from sequence file provided for search

    Arg [3]     : Window size of CRISPR target provided

    Example     : blast(\%CRISPRS, $SEQ, $WINDOWSIZE)

    Description : Run BLAST+ search for CRISPR targets

    Returntype  : Hash reference

    Status      : Development

=cut
sub blast {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '(\%CRISPRfile, $SEQ, $WINDOWSIZE)';
    @_ == 3 or confess wrongNumberArguments(), $filledUsage;

    my ($CRPfile, $seqFile, $word_size) = @_;
    my (%targets, $info);
    $word_size = sprintf "%.0f", ($word_size/2); #round up word_size
    my $BLASTCMD = "blastn -query $CRPfile -subject $seqFile -word_size $word_size -outfmt \"6 qseqid qseqid qstart qend sstart send sstrand pident nident\""; #use 'blastn-short' settings for sequences shorter than 30 nucleotides

    open(BLAST, "$BLASTCMD |") or die "Can't open BLAST commmand <$BLASTCMD>", $!;
    while ( my $blastResult = <BLAST> ) {
        my ($nident) = $blastResult =~ /(\d+)$/; #get number of identical matches
        next if ($nident < $word_size); #skip if match has low identity matches ( < half of $WINDOWSIZE )

        my @result = split('\t', $blastResult);
        my $crispr = $result[0];

        $info = { #anonymous hash with BLAST info for each match
            'qseqid'    => $result[1],
            'qstart'    => $result[2],
            'qend'      => $result[3],
            'sstart'    => $result[4],
            'send'      => $result[5],
            'sstrand'   => $result[6],
            'pident'    => $result[7],
            'nident'    => $result[8],
        };
        # Hash of Array of Hashes to store BLAST results for each query
        # -- Hash key == CRISRP name
        # -- Array accounts for multiple hits for each CRISPR sequence
        # -- Hash contains BLAST match info
        push @{$targets{$crispr}} , $info;
    } close BLAST;

    return(\%targets);
}


#-------------------------------------------------------------------------------
# HELPERS


=head1 COPYRIGHT AND LICENSE

Andres Breton © 2016

[LICENSE]

=head1 CONTACT

Please email comments or questions to Andres Breton, <dev@andresbreton.com>

=head1 SETTING PATH

If PERL5LIB was not set, do something like this:

use FindBin; use lib "$FindBin::RealBin/lib";

This finds and uses subdirectory 'lib' in current directoy as library location

=cut
1;
