package Search;

use Exporter qw(import);
our @ISA = qw(Exporter);
our @EXPORT = qw(findOligo); #functions exported by default
our @EXPORT_OK = qw(blast); #functions for explicit export

use strict; use warnings; use diagnostics; use feature qw(say);
use Carp;

use Bio::Seq; use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlastPlus;
use MyConfig; use MyIO;


# ==============================================================================
#
#   CAPITAN:        Andres Breton, http://andresbreton.com
#   FILE:           Search.pm
#   LICENSE:
#   DEPENDENCIES:   - NCBI's BLAST+ CL utility
#
# ==============================================================================

=head1 NAME

Search - package searching gRNA sequences and offsite targets for CRISPR targeting.

=head1 SYNOPSIS

use Search;

=head1 DESCRIPTION



=head1 EXPORTS

=head2 Default Behaviors

Exports findOligo subroutine by default.

use Search;

=head2 Optional Behaviors

Search::blast;

=head1 FUNCTIONS

=cut

#-------------------------------------------------------------------------------
# MAIN

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
=head2 findOligo

  Arg [1]     : Sequence to be searched details hash

  Arg [2]     : Window size of gRNA sequence

  Example     : findOligo($seqDetails, $windowSize)

  Description : Find gRNA sequences for CRISPR targeting

  Returntype  : Hash of hashes reference

  Status      : Stable

=cut
sub findOligo {
  my $filledUsage = 'Usage: ' . (caller(0))[3] . '($seqDetails, $windowSize)';
  @_ == 2 or confess wrongNumberArguments(), $filledUsage;

  my ($seqDetails, $windowSize) = @_;
  my %gRNAs;
  my $instance = 0; # track gRNA count

  # Anonymous subroutine to find gRNA sequences in forward and reverse strands.
  # Stores both strand findings in same hash (%gRNAs) containing all information
  my $go = sub {
      my ($sequence, $strand) = @_;

      say "Searching gRNA sequences on $strand strand";

      my $seqLen = length($sequence);
      my ($kmer, $PAM, $gRNA, $content, $contentG, $contentC, $GC);

      for (my $i = 0; $i < $seqLen; $i++) {
        my $window = substr $sequence, $i, $windowSize;

        # LAST STEP: When DONE LOOKING UP -- Return gRNA sequences and information.
        # --> Don't go out of bounds when at end of sequence
        # Returns references of HoH containing all gRNA instances found and
        # respective information for each gRNA sequence
        return(\%gRNAs) if ( length($window) < $windowSize );

        if ($window =~ /(.+)(.GG)$/) {
          ($kmer, $PAM) = ($1, $2); # get first 'kmer' and PAM nucleotides in gRNA (kmer + PAM)
          my $name        = "gRNA_$instance"; $instance++;
          my $gRNA        = $kmer . $PAM; # kmer + PAM = gRNA sequence
          my $palindrome  = _palindrome($kmer);

          # GC Content
          $contentG   = $gRNA =~ tr/G//;
          $contentC   = $gRNA =~ tr/C//;
          $GC         = ($contentG + $contentC)/$windowSize;

          # Store gRNA sequences and info in Hash of Hashes
          $content = { # anonymous hash of relevant gRNA content
            'kmer'        => $kmer,
            'palindrome'  => $palindrome,
            'strand'      => $strand,
            'start'       => $i,
            'sequence'    => $gRNA,
            'PAM'         => $PAM,
            'G'           => $contentG,
            'C'           => $contentC,
            'GC'          => $GC,
          };
          # Hash key == gRNA sequence name
          # Hash value == Hash with gRNA content info
          $gRNAs{$name} = $content;
        }
      }
  };

  # Get all gRNA sequences in forward and reverse strands of sequence passed, -seq
  $go->( $seqDetails->{'sequence'}, 'plus' );
  $go->( $seqDetails->{'reverse'}, 'reverse' );
}
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
=head2 blast

  Arg [1]     : gRNA fasta file

  Arg [2]     : Array reference with subject sequence file(s) provided for search

  Arg [3]     : Output file name

  Example     : Search::blast($fastaFile, \@SUBJSEQS, $OUTFILE);

  Description : Run BLAST+ search for gRNA sequences

  Returntype  : Hash reference

  Status      : Stable

=cut
sub blast {
  my $filledUsage = 'Usage: ' . (caller(0))[3] . '(\%fastaFile, \@SUBJSEQS, $OUTFILE)';
  @_ == 3 or confess wrongNumberArguments(), $filledUsage;

  my ($fasta, $SUBJSEQS, $OUTFILE) = @_;
  my (%targets, $info, $hsps);
  my @SUBJSEQS  = @$SUBJSEQS;
  my $outDir    = $main::OUTDIR;
  my $wordSize  = 7;

  mkDir("$outDir/blast");

  foreach my $subject (@SUBJSEQS) {
    my $subjName  = _getSeqName($subject);
    my $outFile   = "$outDir/blast/$subjName\_$OUTFILE\_blast.txt";

    say "Searching gRNA sequences against $subject";

    # Create StandAloneBlastPlus Factory
    my $fac = Bio::Tools::Run::StandAloneBlastPlus->new(
          -db_data  => $subject,
        );

    # Perform BLAST call
    $fac->blastn( -query => $fasta, -outfile => $outFile,
                  -method_args => [ -word_size => 7],
          );

    # Rewind to beginning of results and get all
    $fac->rewind_results;

    # Process each gRNA sequence
    while ( my $result = $fac->next_result ) {
      my ($gRNA)  = $result->query_name =~ /(.*):\d+/;; # gRNA sequence name ex.) 'gRNA_0', removes appendend positioning
      my $numHits = $result->num_hits;

      # Resolve when gRNA sequences has no matches in subject
      say "\t$gRNA has $numHits hits here" if ($numHits == 0);

      # Default anonymous hashes for storing BLAST results
      # $info hash stores each gRNA BLAST result's info
      # $hsps hash stores that match's High-scoring Segment Pairs (HSPs) result
      $info = {
        'numhits'     => $numHits,
        'occurrences' => 0,
      };

      $hsps = {
        'rank'    => 0,
        'qstart'  => 0,
        'qend'    => 0,
        'sstart'  => 0,
        'send'    => 0,
        'sstrand' => 0,
        'pident'  => 0,
        'nident'  => 0,
        'gaps'    => 0,
      };

      # Process each gRNA hit
      while ( my $hit = $result->next_hit ) {

        $info->{'occurrences'} = $hit->num_hsps;

        # Process each match (HSP) in iterative fashion
        while( my $hsp = $hit->next_hsp ) {
          # Get all values to store in hashes
          $hsps = {
            'rank'    => $hsp->rank,
            'qstart'  => $hsp->start('query'),
            'qend'    => $hsp->end('query'),
            'qstrand' => $hsp->strand('query'),
            'sstart'  => $hsp->start('hit'),
            'send'    => $hsp->end('hit'),
            'sstrand' => $hsp->strand('hit'),
            'pident'  => $hsp->percent_identity,
            'nident'  => $hsp->num_identical,
            'gaps'    => $hsp->gaps,
          };
          push @{ $targets{$gRNA}{$subjName}{'hsps'} } , $hsps;
        }
      }
      push @{ $targets{$gRNA}{$subjName}{'info'} } , $info;
    }
    $fac->cleanup; #clean up temp database files
  }

  say "\nBLAST files saved in: '$outDir/blast' ";

  # Hash of Hashes of Hashes of Arrays of Hash to store BLAST results for each query
  # -- Hash key == gRNA name
  # -- Hash key == Subject name
  # -- Hash key == 'info'
  # -- Hash key == 'hsps'
  # -- Array accounts for multiple hits for each gRNA sequence as hashes....
  # -- Hash contains BLAST match info
  return(\%targets);
}

#-------------------------------------------------------------------------------
# HELPERS

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ($seq);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes 1 argument, a record file and extracts relevant information
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = ($name); Name of sequence record
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub _getSeqName {
  my ($seq) = @_;

  # Sequence OO
  my $seqInObj  = Bio::SeqIO->new(-file => $seq, -alphabet => "dna");
  my $format    = $seqInObj->_guess_format($seq); #check format of input file
  my $seqObj    = $seqInObj->next_seq;
  my $name      = $seqObj->display_id;

  return $name;
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ($sequence);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes 1 argument, a sequence string and determines if string is
# palindrome.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = Returns yes/no if palindrome found
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

sub _palindrome {
  my $filledUsage = 'Usage: ' . (caller(0))[3] . '($sequence)';
  @_ == 1 or confess wrongNumberArguments(), $filledUsage;

  my ($sequence) = @_;

  my $palindrome = reverse $sequence;

  $palindrome eq $sequence ? return("Yes") : return("No");
}



=head1 COPYRIGHT AND LICENSE

Andres Breton (C) 2017

[LICENSE]

=head1 CONTACT

Please email comments or questions to Andres Breton, dev@andresbreton.com

=head1 SETTING PATH

If PERL5LIB was not set, do something like this:

use FindBin; use lib "$FindBin::RealBin/lib";

This finds and uses subdirectory 'lib' in current directory as library location

=cut
1;
