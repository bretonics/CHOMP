package Parser;

use Exporter qw(import);
our @ISA = qw(Exporter);
our @EXPORT = qw(parseHeader parseFeatures); #functions exported by default

use warnings; use strict; use diagnostics; use feature qw(say);
use Carp;
use Bio::DB::GenBank; use Bio::SeqFeatureI;

# =============================================
#
#   CAPITAN: Andres Breton
#   FILE: parser.pm
#
# =============================================

#-------------------------------------------------------------------------
# MAIN

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ();
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes  arguments
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub parseHeader {
    my ($NCBIfile) = @_;
    unless(open(INFILE, "<", $NCBIfile)) {
        croak "Can't open $NCBIfile for reading " , $!;
    }
    my ($locus, $seqLen, $accession, $version, $gi, $organism, $sequence) = qw(NA NA NA NA NA NA NA);

    # Slurp File
    $/ = ''; #line separator
    my $FILE = <INFILE>;
    $/ = "\n";  #set back line separator
    close INFILE;

    say "Getting [header] content...";

    # Get Locus Name and Sequence Length
    ($locus, $seqLen) = getLocus($FILE);
    # Get Accession
    $accession = getAccession($FILE);
    # Get Version
    $version = getVersion($FILE);
    # Get GI
    $gi = getGI($FILE);
    # Get Organism
    $organism = getOrganism($FILE);

    return $locus, $seqLen, $accession, $version, $gi, $organism;
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ();
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes  arguments
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub parseFeatures {
    my ($id) = @_;
    my ($proteinID, $translation, $gene, $exon) = qw(NA NA NA NA);
    my $dbObject = Bio::DB::GenBank->new;   #set database object
    my $seqObject = $dbObject->get_Seq_by_id($id);  #set seq object
    say "Getting [features] content...\n";
    my $sequence = $seqObject->seq();
    for my $feature ($seqObject->get_SeqFeatures) {   #gets seqObject features
        # Get Protein ID and Translation
        if($feature->primary_tag eq "CDS") {
            ($proteinID, $translation) = getProteinID($feature);
        }
        # Get Exon
        if ($feature->primary_tag eq "exon") {
            ($exon) = getExon($feature);
        }
        # Get Gene
        if ($feature->primary_tag eq "gene") {
            ($gene) = getGene($feature);
        }
    }
    return $sequence, $proteinID, $translation, $gene;
}

#-------------------------------------------------------------------------
# HELPERS
sub getLocus {
    my ($file) = @_;
    if($file =~ /^LOCUS\s+(\w+)\s+(\d+)\s+/) {
        my ($locus, $seqLen) = ($1, $2);
        return $locus, $seqLen;
    }
}

sub getAccession {
    my ($file) = @_;
    if($file =~ /^ACCESSION\s+(\w+)/m) {
        my $accession = $1; return $accession;
    }
}

sub getVersion {
    my ($file) = @_;
    if($file =~ /^VERSION\s+(\w+)\s+/m) {
        my $version = $1; return $version;
        }
}

sub getGI {
	my ($file) = @_;
	if($file =~ /^VERSION.*GI:(\w+)/m){
		my $gi = $1; return $gi;
	}
	else{
		croak "ERROR getting GI", $!;
	}
}

sub getOrganism {
    my ($file) = @_;
    if($file =~ /organism="(.*?)"/) {
        my $organism = $1; return $organism;
    }
}

sub getSequence {
    my ($file) = @_;
    my $seq;
	if($file =~ /ORIGIN\s*(.*)\/\//s){
	    $seq = $1;
	}
	else{
		croak "ERROR getting sequence";
	}
	$seq =~ s/[\s\d]//g; #remove spaces/numbers from sequence
	return uc($seq);
}

sub getGene {
	my ($feature) = @_;
    my $gene;
	if($feature->has_tag("gene")) {
        ($gene) = $feature->get_tag_values("gene");
	}else{
		# return "unknown";
	}
    return $gene;
}

sub getExon {
    my ($feature) = @_;
    my $exon;
    if($feature->has_tag("exon")) {
        ($exon) = $feature->get_tag_values("exon");
	}else{
		# return "unknown";
	}
    return $exon;
}

sub getProteinID {
	my ($feature) = @_;
    my ($proteinID, $translation);
    if($feature->has_tag("protein_id")) {
        ($proteinID) = $feature->get_tag_values("protein_id");
        ($translation) = $feature->get_tag_values("translation");
    }else{
		# return "unknown";
	}
    return $proteinID, $translation;
}
1;
