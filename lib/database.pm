package Database;

use Exporter qw(import);
our @ISA = qw(Exporter);
our @EXPORT = qw(startMongoDB insertData updateData readData removeData); #functions exported by default

use warnings; use strict; use diagnostics; use feature qw(say);
use Carp;

use MongoDB; use MongoDB::OID;


# =============================================
#
#   CAPITAN: Andres Breton
#   FILE: database.pm
#
# =============================================

my @dataFields = qw(_id accession sequence version locus organism seqLength gene proteinID translation);
#-------------------------------------------------------------------------
# MAIN
sub startMongoDB {
    my ($MONGODB, $outDir) = @_;
    my $pid;
    my $dbDir = $outDir."/db";
    `mkdir $dbDir` unless (-e $dbDir);
    my $mongoLog = $dbDir."/mongo.log";

    my $command = "mongod --dbpath $dbDir --logpath $mongoLog --fork";
    say "\nStarting MongoDB server...";
    my @result = `$command`; #get shell results
    if ($? == 0) { #Check return value
        $pid = $result[1] =~ /.+:\s(\d+)$/; $pid = $1; #get child PID
        say "MongoDB successfully started.\n";
        return $pid;
    } elsif ($? == 25600) { #Possible mongd already running
        say "*********FAILED";
        say "Could not fork. This was most likely caused by an instance of [mongod] already running.";
        # Check for Currently Running MongoDB Server
        my @mongdPS = `ps -e -o pid,args | grep \"mongod\"`;
        if ($mongdPS[0] =~ /^\s?(\d+)\s+mongod.*/) {
            $pid = $1;
            say "YES! Found running process: $mongdPS[0]";
            print "Would you like to continue (y/n)? ";
            my $response = lc <>; chomp $response;
            if ($response eq "yes" || $response eq "y") {
                return $pid;
            } else {
                exit;
            }
        } else {
            croak "Sorry, could not find instance of mongod running on system. Please check processes.", $!;
        }
    } else {
        croak "ERROR: Failed to execute $command\n Something happened that did not allow MongoDB server to start!", $!;
    }
}

sub insertData {
    my ($MONGODB, $COLLECTION, $id, $gi, $accession, $version, $locus, $organism, $sequence, $seqLen, $gene, $proteinID, $translation) = @_;

    my $collectionObj = databaseConnection($MONGODB, $COLLECTION);
    say "Storing data for ID ($id) into database $MONGODB";
    $collectionObj->insert({_id => $gi, #GI stored as Mongo UID
                        "accession" => $accession,
                        "version" => $version,
                        "locus" => $locus,
                        "organism" => $organism,
                        "seqLength" => $seqLen,
                        "sequence" => $sequence,
                        "gene" => $gene,
                        "proteinID" => $proteinID,
                        "translation" => $translation
                        })
}

sub updateData {
    my ($field, $value, $MONGODB, $COLLECTION) = @_;
    say "\nUPDATING $field record [$value] in database...";
    say "Available fields are:\t@dataFields\n";
    print "What field do you want? ";
    my $fieldUpdate = <>; chomp $fieldUpdate;
    print "What is the NEW value for $fieldUpdate field? ";
    my $newValue = <>; chomp $newValue;
    my $collectionObj = databaseConnection($MONGODB, $COLLECTION);
    $collectionObj->update({$field => $value}, {'$set' => {$fieldUpdate => $newValue}});
    say "Document $value updated. $fieldUpdate field changed to $newValue.";
}

sub readData {
    my ($field, $value, $MONGODB, $COLLECTION) = @_;
    my $collectionObj = databaseConnection($MONGODB, $COLLECTION);
    say "\nREADING field \"$field\" value \"$value\" from database...";
    my $cursor = $collectionObj->find({$field => $value});
    while (my $obj = $cursor->next) {
        say "Available fields are:\t@dataFields\n";
        print "What field do you want? ";
        my $response = <>; chomp $response;
        say "Here you go [$response]:\n", $obj->{$response};
    }
}

sub removeData {
    my ($field, $value, $MONGODB, $COLLECTION) = @_;
    say "REMOVING $field record [$value] in database...";
    my $collectionObj = databaseConnection($MONGODB, $COLLECTION);
    $collectionObj->remove({$field => $value});
}
-------------------------------------------------------------------------
# HELPERS
sub databaseConnection {
    my ($MONGODB, $COLLECTION) = @_;
    my $client = MongoDB::MongoClient->new; #connect to local db server
    my $db = $client->get_database($MONGODB); #get MongoDB databse
    my $collectionObj = $db->get_collection($COLLECTION); #get collection
    return $collectionObj;
}
1;
