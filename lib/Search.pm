package Search;

use Exporter qw(import);
our @ISA = qw(Exporter);
our @EXPORT = qw(blast); #functions exported by default
our @EXPORT_OK = qw(); #functions for explicit export

use strict; use warnings; use diagnostics; use feature qw(say);
use Carp;

use MyConfig;

# ==============================================================================
#
#   CAPITAN: Andres Breton http://andresbreton.com
#   FILE: search.pm
#   LICENSE:
#   USAGE:
#   DEPENDENCIES:   - NCBI's BLAST+ CL utility
#
# ==============================================================================

=head1 NAME

Search - package calling a BLAST search to find CRISPR offsite targets

=head1 SYNOPSIS

Creation:
    use Search;

=head1 DESCRIPTION


=head1 EXPORTS

=head2 Default Behaviors

Exports $SUB subroutine by default

use Search;

=head2 Optional Behaviors

Search::;

=head1 FUNCTIONS

=cut

#-------------------------------------------------------------------------------
# MAIN

=head2 blast

    Arg [1]     :

    Example     :

    Description : Run BLAST+ search for CRISPR targets

    Returntype  :

    Status      : Development

=cut
sub blast {
    my $filledUsage = 'Usage: ' . (caller(0))[3] . '(\@CRISPRS)';
    @_ == 1 or confess wrongNumberArguments(), $filledUsage;

    my (@CRISPRS) = @_;

    return;
}

#-------------------------------------------------------------------------------
# HELPERS


=head1 COPYRIGHT AND LICENSE

Andres Breton © 2016

[LICENSE]

=head1 CONTACT

Please email comments or questions to Andres Breton me@andresbreton.com

=head1 SETTING PATH

If PERL5LIB was not set, do something like this:

use FindBin; use lib "$FindBin::RealBin/lib";

This finds and uses subdirectory 'lib' in current directoy as library location

=cut
1;
