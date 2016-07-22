#!/usr/bin/perl

use strict; use warnings; use diagnostics; use feature qw(say);
use Getopt::Long; use Pod::Usage;

use File::Basename qw(dirname);
use Cwd qw(abs_path);
use FindBin; use lib "$FindBin::RealBin/lib";

use Eutil; use Parser; use Database;

use Bio::Seq; use Bio::SeqIO;

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
#
#   CAPITAN: Andres Breton
#   FILE: chomp.pl 🐊
#   LICENSE:
#   USAGE:
#   DEPENDENCIES:
#
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# USER VARIABLES
my $SP6 = "ATTTAGGTGACACTATA";
my $OVERLAP = "GTTTTAGAGCTAGAAATAGCAAG";
my $CONSTOLIGO = "AAAAGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC";
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# COMMAND LINE
my $FILE = "";
my @IDS;
my $MONGODB = "CHOMP"; #defaults to
my $COLLECTION = "";
my ($INSERT, @UPDATE, @READ, @REMOVE);
my $USAGE= "\n\n $0 [options]\n
Options:
    -file           File
    -ids            ID(s)
    -mongo          MongoDB database name
    -collection     Collection name in MongoDB database
    -insert         Insert into database [optional/default]
    -update         Update database
    -read           Read from database
    -remove         Remove from database
    -help           Shows this message
\n";

# OPTIONS
GetOptions(
    'file:s'        =>\$FILE,
    'id:i{1,}'      =>\@IDS,
    'mongo:s'       =>\$MONGODB,
    'mongo:s'       =>\$MONGODB,
    'collection:s'  =>\$COLLECTION,
    'insert+'       =>\$INSERT,
    'update:s{1,}'  =>\@UPDATE,
    'read:s{1,}'    =>\@READ,
    'remove:s{1,}'  =>\@REMOVE,
    help            =>sub{pod2usage($USAGE);}
)or pod2usage(2);
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# VARIABLES
my $REALBIN = "$FindBin::RealBin";
my $OUTDIR = createOutDir("Data");

#Color Output
my $GRNTXT = "\e[1;32m"; #bold green
my $REDTXT = "\e[1;31m"; #bold red
my $NC = "\e[0m"; #color reset


my $seqInObject = Bio::SeqIO->new(-file => "$inFile", -format => "genbank", -alphabet => "dna");
my $format = $seqInObject->_guess_format("$inFile"); #check format of input file

my $sequence = $seqInObject->next_seq;
my $actual = $sequence->seq;

my ($fileName) = $FILE =~ /(\w+)\b\./; #extract file name
my $seqOutObject = Bio::SeqIO->new(-file => ">$fileName.fasta", -format => "fasta", -alphabet => "dna");
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CALLS
checks();

#lookUpFeatures($seqObject);

$seqOutObject->write_seq($sequence);    #write FASTA file from input file

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SUBS

#---------------------------------------------------------
# $input = checks();
#---------------------------------------------------------
# This function checks for arguments passed in command-line
# using global variables
#---------------------------------------------------------
#
#---------------------------------------------------------
sub checks {
    unless ($FILE){
        die "Did not provide an input file, -file <infile.txt>", $USAGE;
    }
    if ($format ne "genbank") {
        say "Incorrect file format. You provided $format format.";
        say "Please provide a GenBank file format";
    }
    my $argLen = @ARGV;
    unless ($argLen >= 1) {
        say "Did not provide the correct number of paramaters", $USAGE;
    }
}

#---------------------------------------------------------
# $input = ();
#---------------------------------------------------------
# This function takes  arguments
#
#---------------------------------------------------------
#
#---------------------------------------------------------
sub lookUpFeatures {
    my (@seqObjects) = @_;
    my $arraySize = @seqObjects;
    for(my $i=0; $i<$arraySize;$i++) {  #loop through seqObjects passed
        for my $feat ($seqObjects[$i]->get_SeqFeatures) {   #gets seqObject features
            # Get Protein ID and Translation
            if ($feat->primary_tag eq "CDS") {
                getFeatures($feat, $feat->primary_tag);
            }
            # Get Exon
            if ($feat->primary_tag eq "exon") {
                getFeatures($feat, $feat->primary_tag);
            }
            }

        }
    }
}

#---------------------------------------------------------
# $input = ();
#---------------------------------------------------------
# This function takes  arguments
#
#---------------------------------------------------------
#
#---------------------------------------------------------
sub getFeatures {
    my ($feat, $primaryTag) = @_;
    print "\nPrimary Tag: ", $feat->primary_tag, " start: ", $feat->start, " ends: ", $feat->end, " strand: ", $feat->strand,"\n";
    for my $tag ($feat->get_all_tags) { #gets seqObject tags from primary feature
        print " tag: ", $tag, "\n";
        for my $value ($feat->get_tag_values($tag)) { #gets seq object values from tag
            print "  value: ", $value, "\n";
        }
    }

}
#---------------------------------------------------------
# $result = ();
#---------------------------------------------------------
# This function takes  arguments
#
#---------------------------------------------------------
#
#---------------------------------------------------------
sub createOutDir {
    my $outDir = @_;
    `mkdir $outDir` unless(-e $outDir);
    return $outDir;
}
sub createDownloadDir {
    my $outDir = "Data";
    if (! -e $outDir){
        `mkdir $outDir`;
    }
    return $outDir; subhd
}

sub compare {
    my () = @_;

}
