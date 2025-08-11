package SS;

use Exporter qw(import);
our @ISA = qw(Exporter);
our @EXPORT = qw(ss); # functions exported by default
our @EXPORT_OK = qw(); # functions for explicit export

use strict; use warnings; use diagnostics; use feature qw(say);
use Carp;

use MyConfig; use MyIO;

use RNA; # ViennaRNA

# ==============================================================================
#
#   CAPITAN:        Andres Breton
#   FILE:           SSpm
#   LICENSE:        [LICENSE]
#   USAGE:          Predict RNA secondary structures
#   DEPENDENCIES:   ViennaRNA (https://www.tbi.univie.ac.at/RNA/)
#
# ==============================================================================

=head1 NAME

SS - package for prediction and comparison of RNA secondary structures using ViennaRNA.

=head1 SYNOPSIS

use SS;

=head1 DESCRIPTION

This module is designed to predictic secondary structures of nucleotide sequences using ViennaRNA package.

=head1 EXPORTS

=head2 Default Behaviors

Exports SS subroutine by default.

use SS;

=head1 FUNCTIONS

=cut

#-------------------------------------------------------------------------------
# MAIN

=head2 ss

  Arg [1]     : Name string

  Arg [2]     : Sequence string

  Example     : ss($name, $sequence)

  Description : Predict secondary structures and output plots

  Returntype  : undef

  Status      : Stable

=cut
sub ss {
  my $filledUsage = 'Usage: ' . (caller(0))[3] . '($name, $sequence)';
  @_ == 2 or confess wrongNumberArguments(), $filledUsage;

  my ($name, $sequence) = @_;
  my $outFile = "$main::OUTDIR/ss/$name";
  my $rnaPlot = "$outFile\_rna.ps";
  my $dotPlot = "$outFile\_dot.ps";
  my $rssPlot = "$outFile\_rss.ps";

  my $F = RNA::pf_fold($sequence);   # compute partition function and pair pobabilities
  my ($ss, $mfe) = RNA::fold($sequence);
  RNA::PS_rna_plot($sequence, $ss, $rnaPlot);  # write PS plot to gRNA_#_rna.ps
  RNA::PS_dot_plot($sequence, $dotPlot);     # write dot plot to dot.ps

  my $command = "relplot.pl $rnaPlot $dotPlot > $rssPlot";
  say "Running relative secondary structure command --> '$command'";
  `$command`;

  return;
}

=head1 COPYRIGHT AND LICENSE

Andres Breton (C)

[LICENSE]

=head1 CONTACT

Please email comments or questions to Andres Breton, dev@andresbreton.com

=head1 SETTING PATH

If PERL5LIB was not set, do something like this:

use FindBin; use lib "$FindBin::RealBin/lib";

This finds and uses subdirectory 'lib' in current directory as library location

=cut

1;
