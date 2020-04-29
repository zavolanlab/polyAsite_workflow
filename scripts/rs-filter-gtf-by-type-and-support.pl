use strict;
use warnings;
use Getopt::Long;

my @types = ();
my $type_id;
my $support_level;
my $support_level_id;

GetOptions(
  "type=s"             => \@types,
  "type_id=s"          => \$type_id,
  "support_level=i"    => \$support_level,
  "support_level_id=s" => \$support_level_id
);

if ( $#types == -1 or not defined $ARGV[0] ) {
  print STDERR
    "[ERROR] Usage: perl $0 --type_id=transcript_biotype --type=protein_coding --type=lincRNA Homo_sapiens.GRCh38.85.gtf.gz\n";
  print STDERR
    "[ERROR] If support level is given via --support_level, only transcript with this support level or lower are considered\n";
  print STDERR "[ERROR] In this case, also --support_level_id is required (e.g. 'transcript_support_level')\n\n";
  exit(2);
}

if( (defined $support_level_id && not defined $support_level) ||
    (defined $support_level && not defined $support_level_id) ) {
  print STDERR "[ERROR] Both --support_level and --support_level_id are required if either of both is given\n\n";
  exit(2);
}

my %valid = map { $_ => 1 } @types;

my %genes     = ();
my %processed = ();

open( F, "< $ARGV[0] " ) or die "Can't open $ARGV[0]\n";
while (<F>) {
  if ( $_ =~ /^#/ ) {
    print $_;
    next;
  }

  my @F = split("\t");
  my $id;
  if ( $_ =~ /gene_id\s"([^"]+)"/ ) {
    $id = $1;
  } else {
    print STDERR "[ERROR] Could not infer gene id for:\n$_";
    exit(2);
  }
  if ( "gene" eq $F[2] ) {
    $genes{$id} = $_;
    next;
  }
  my $type;
  if ( $_ =~ /$type_id\s"([^"]+)"/ ) {
    $type = $1;
  } else {
    print STDERR "[ERROR] Could not infer transcript type for:\n$_";
    exit(2);
  }
  my $support;
  if ( defined $support_level && $_ =~ /$support_level_id\s"([^"|\s]+)/ ) {
    $support = $1;
  }

  if ( defined $support_level ) {
    next if ( $support eq "NA" or $support > $support_level );
  }

  if ( defined $valid{$type} ) {
    if ( not defined $processed{$id} ) {
      print $genes{$id};
      $processed{$id} = 1;
    }
    print $_;
  }
}
close(F);
