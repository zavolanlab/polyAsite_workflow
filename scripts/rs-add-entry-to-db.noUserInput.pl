use strict;
use warnings;
use MongoDB;
use MongoDB::OID;
use Getopt::Long;
use boolean;

my $help;
my $SERIESID;
my $SAMPLEID;
my $sample_db_id;
my $adapter5p   = '';
my $adapter3p   = '';
my $prot        = '';
my $sourceInput = "NA";
my $sexInput    = "NA";
my $organism;
my $title;
my $genome_version;
my $reversedInput;
my $treatment = 'none';
my $host;
my $db_name;
my $collection_name;
my $pub_id;

GetOptions(
  "series=s"           => \$SERIESID,
  "sample=s"           => \$SAMPLEID,
  "sample_genome_id=s" => \$sample_db_id,
  "help"               => \$help,
  "h"                  => \$help,
  "5pAd=s"             => \$adapter5p,
  "3pAd=s"             => \$adapter3p,
  "protocol=s"         => \$prot,
  "source=s"           => \$sourceInput,
  "sex=s"              => \$sexInput,
  "treatment=s"        => \$treatment,
  "toReverse=s"        => \$reversedInput,
  "organism=s"         => \$organism,
  "genome=s"           => \$genome_version,
  "host=s"             => \$host,
  "db=s"               => \$db_name,
  "collection=s"       => \$collection_name,
  "pubmed=s"           => \$pub_id,
  "title=s"            => \$title
);

my $showhelp = 0;
$showhelp = 1 if ( not defined $SERIESID );
$showhelp = 1 if ( not defined $SAMPLEID );
$showhelp = 1 if ( defined $help );

if ($showhelp) {
  print STDERR "Usage: $0 --series=GSE27452 --sample=GSM678411\n\n";
  exit;
}

if ( not defined $reversedInput or $reversedInput !~ /^([Yy]|[Nn])/ ) {
  print STDERR "[ERROR] toReverse not set\n";
  print STDERR "[ERROR] Valid strings: Y/N\n";
  exit;
}

if ( $reversedInput =~ /^Y/i ) {
  $reversedInput = true;
} else {
  $reversedInput = false;
}

# Connection
my $connection = MongoDB::Connection->new( host => "mongodb://$host" );
my $user       = '';
my $pw         = '';
$connection->authenticate( $db_name, $user, $pw );

# Database
my $db = $connection->get_database($db_name);

# Collection
my $samples = $db->get_collection($collection_name);

# check if an entry of this sample for this genome version already exists
if ( $samples->find_one( { "sample_genome_id" => $sample_db_id } ) ) {
  print STDERR "[INFO] Sample $SAMPLEID for genome $genome_version already found in DB. Nothing done!\n";
  exit(0);
}

$samples->insert( {
    "GEOseriesID"      => $SERIESID,
    "GEOsampleID"      => $SAMPLEID,
    "sample_genome_id" => $sample_db_id,
    "genome"           => $genome_version,
    "PubMedID"         => $pub_id,
    "title"            => $title,
    "organism"         => $organism,
    "sex"              => $sexInput,
    "reads"            => {
      "raw"     => { "nr" => 0, "md5sum" => '', "length" => { "max" => 0, "min" => 0 } },
      "trimmed" => {
        "nr"              => 0,
        "md5sum"          => '',
        "md5sumcollapsed" => '',
        "md5sumNoPolyA"   => '',
        "md5sum5p"        => "",
        "length"          => { "max" => 0, "min" => 0 }
      },
      "genome" => {
        "uniqueMappers" => {
          "md5sum" => "",
          "nr"     => 0
        },
        "multiMappers" => {
          "md5sum" => "",
          "nr"     => 0
        }
      },
      "mapped" => {
        "uniqueMappers" => {
          "md5sum" => "",
          "nr"     => 0
        }
      },
      "transcriptome" => {
        "uniqueMappers" => {
          "md5sum" => "",
          "nr"     => 0
        }
      },
      "valid" => {
        "nr"      => 0,
        "md5sum"  => '',
        "cutoffs" => {
          "maxLength" => 0,
          "maxAs"     => 0.8,
          "maxNs"     => 2
        }
      },
    },
    "sites" => {
      "all" => {
        "internalpriming" => {
          "md5sum" => "",
          "number" => { "plus" => 0, "minus" => 0 },
          "reads"  => { "plus" => 0, "minus" => 0 }
        },
        "md5sum" => '',
        "number" => { "plus" => 0, "minus" => 0 },
        "reads"  => { "plus" => 0, "minus" => 0 }
      },
      "highconfidence" => {
        "internalpriming" => {
          "md5sum" => "",
          "number" => { "plus" => 0, "minus" => 0 },
          "reads"  => { "plus" => 0, "minus" => 0 }
        },
        "md5sum" => '',
        "number" => { "plus" => 0, "minus" => 0 },
        "reads"  => { "plus" => 0, "minus" => 0 }
      }
    },
    "adaptor3p"    => $adapter3p,
    "adaptor5p"    => $adapter5p,
    "protocol"     => $prot,
    "log"          => [],
    "tobereversed" => $reversedInput,
    "source"       => $sourceInput,
    "treatment"    => $treatment,
    "visible"      => true,
    "pairedEnd"    => false
  }
);

