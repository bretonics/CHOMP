package SS;

use Exporter qw(import);
our @ISA = qw(Exporter);
our @EXPORT = qw(ss); # functions exported by default
our @EXPORT_OK = qw(); # functions for explicit export

use strict; use warnings; use diagnostics; use feature qw(say);
use Carp;

use MyConfig; use MyIO;

use RNA; #ViennaRNA

sub ss {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '($name, $sequence)';
    @_ == 2 or confess wrongNumberArguments(), $filledUsage;

    my ($name, $sequence) = @_;
    my $outFile = "ss/$name";
    my $rnaPlot = "$outFile\_rna.ps";
    my $dotPlot = "$outFile\_dot.ps";
    my $rssPlot = "$outFile\_rss.ps";
    my $F = RNA::pf_fold($sequence);   # compute partition function and pair pobabilities
    my ($ss, $mfe) = RNA::fold($sequence);
    RNA::PS_rna_plot($sequence, $ss, $rnaPlot);  # write PS plot to CRISPR_#_rna.ps
    RNA::PS_dot_plot($sequence, $dotPlot);       # write dot plot to dot.ps

    my $command = "relplot.pl $rnaPlot $dotPlot > $rssPlot";
    say "Running relative secondary structure command '$command'";
    `$command`;

    return;
}



1;
