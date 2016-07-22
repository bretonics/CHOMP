package Eutil;

use Exporter qw(import);
our @ISA = qw(Exporter);
our @EXPORT = qw(getNCBIfile); #functions exported by default

use warnings; use strict; use diagnostics; use feature qw(say);
use Carp;
use Bio::DB::EUtilities;

# =============================================
#
#   CAPITAN: Andres Breton
#   FILE: eutil.pm
#
# =============================================

#-------------------------------------------------------------------------
# MAIN
sub getNCBIfile {
    my ($ID, $outDir, $FORCE, $DATABASE, $TYPE, $email) = @_;
    my $outFile = $outDir ."/$ID." .$TYPE;

    my $eutil = Bio::DB::EUtilities->new(
            -eutil => "efetch",
            -db => $DATABASE,
            -id => $ID,
            -email => $email,
            -rettype => $TYPE,
    );

    #Fetch
    say "\nFetching data from NCBI for ID $ID. Wait...\n";
    $eutil->get_Response( -file => $outFile);
    sleep(3); #Don't overload NCBI requests
    return ($outFile, 1);
}
1;
