package Search;

use Exporter qw(import);
our @ISA = qw(Exporter);
our @EXPORT = qw(findOligo); #functions exported by default
our @EXPORT_OK = qw(blast); #functions for explicit export

use strict; use warnings; use diagnostics; use feature qw(say);
use Carp;

use MyConfig; use MyIO;

use Data::Dumper;

# ==============================================================================
#
#   CAPITAN:        Andres Breton, http://andresbreton.com
#   FILE:           Search.pm
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

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
=head2 findOligo

    Arg [1]     : Sequences info hash to be searched

    Arg [2]     : Window size of CRISPR target

    Example     : findOligo($seqInfo, $windowSize)

    Description : Find CRISPR targets

    Returntype  : Hash of hashes reference

    Status      : Development

=cut
sub findOligo {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '($seqInfo, $windowSize)';
    @_ == 2 or confess wrongNumberArguments(), $filledUsage;

    my ($seqInfo, $windowSize) = @_;
    my (%CRISPRS, @CRPseqs);
    my $instance = 0; # track CRISPR count

    # Anonymous subroutine to find CRISPR sequences in forward and reverse strands.
    # Stores both strand findings in same hash (%CRISPRS) containing all information
    # and array (@CRPseqs) containing all CRISPR sequences [cause why not].
    my $go = sub {
            my ($sequence, $strand) = @_;

            my $seqLen = length($sequence);
            my ($gRNA, $PAM, $content, $contentG, $contentC, $GC);

            for (my $i = 0; $i < $seqLen; $i++) {
                my $window = substr $sequence, $i, $windowSize;

                # LAST STEP: When DONE LOOKING UP -- Return CRISPR sequences and information
                if ( length($window) < $windowSize ) { #don't go out of bounds when at end of sequence
                    foreach my $name (keys %CRISPRS) {
                        my $crispr = $CRISPRS{$name}{'gRNA'} . $CRISPRS{$name}{'PAM'}; #join gRNA + PAM sequence
                        $CRISPRS{$name}{'sequence'} = $crispr; # add CRISPR sequence (gRNA + PAM) to each hash
                        push @CRPseqs, $crispr # push each N-oligomer CRISPR seq to array
                    }
                    # Return references of HoH containing all CRISPR instances found and respective information for each + array with the all CRISPR sequences joined (kmer gRNA + PAM)
                    return(\%CRISPRS, \@CRPseqs);
                };

                if ($window =~ /(.+)(.GG)$/) {
                    ($gRNA, $PAM) = ($1, $2); # get first 'kmer' number of nucleotides in gRNA (kmer) + PAM (NGG), gRNA + PAM = crispr sequence
                    my $name = "CRISPR_$instance"; $instance++;

                    # GC Content
                    $contentG   = $window =~ tr/G//;
                    $contentC   = $window =~ tr/C//;
                    $GC         = ($contentG + $contentC)/$windowSize;

                    # Store CRISPR oligomers and info in Hash of Hashes
                    $content = { #anonymous hash of relevant gRNA content
                        'strand'    => $strand,
                        'start'     => $i,
                        'gRNA'      => $gRNA,
                        'PAM'       => $PAM,
                        'G'         => $contentG,
                        'C'         => $contentC,
                        'GC'        => $GC,
                    };
                    # Hash key == CRISPR sequence name
                    # Hash value == Hash with CRISPR content info
                    $CRISPRS{$name} = $content;
                }
            }
    };

    # Get all CRISPR sequences in forward and reverse strands of sequence passed, -seq
    $go->( $seqInfo->{'sequence'}, 'plus' );
    $go->( $seqInfo->{'reverse'}, 'reverse' );
}
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '(\%CRISPRfile, $SUBJSEQ, $WINDOWSIZE, $HTML)';
    @_ == 4 or confess wrongNumberArguments(), $filledUsage;

    my ($CRPfile, $SUBJSEQ, $wordSize, $HTML) = @_;
    my (%targets, $info);
    $wordSize = 7;
    # $wordSize = sprintf "%.0f", ($wordSize/2); #Make wordSize == 1/2 of WINDOWSIZE when searching BLAST hits

    # Use 'blastn-short' settings for sequences shorter than 30 nucleotides
    my $BLASTCMD = "blastn -query $CRPfile -subject $SUBJSEQ -word_size $wordSize -outfmt \"6 qseqid sseqid qstart qend sstart send sstrand pident nident\"";

    open(BLAST, "$BLASTCMD |") or die "Can't open BLAST commmand <$BLASTCMD>", $!;
    while ( my $blastResult = <BLAST> ) {
        my ($nident) = $blastResult =~ /(\d+)$/; # get number of identical matches
        # next if ($nident < $wordSize); #skip if match has low identity matches ( < half of $WINDOWSIZE )

        my @result = split('\t', $blastResult);
        my ($crispr) = $result[0] =~ m/(.+):\d+/; # CRISPR sequence name ex.) 'CRISPR_0' instead of 'CRISPR_0:0'

        $info = { #anonymous hash with BLAST info for each match
            'sseqid'    => $result[1],
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
        push @{ $targets{$crispr} } , $info;
    } close BLAST;

    my $BLASTCMD_HTML = "blastn -query $CRPfile -subject $SUBJSEQ -word_size $wordSize -out blast.html -html";
    `$BLASTCMD_HTML` if($HTML);

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
