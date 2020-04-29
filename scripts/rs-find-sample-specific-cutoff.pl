use strict;
use warnings;
use Getopt::Long;
#use Memory::Usage;

#my $mu = Memory::Usage->new();
#$mu->record('starting work');

my $cutoff;
my $sample;
my $upstream;
my $downstream;

GetOptions(
  "cutoff=i"     => \$cutoff,
  "sample=s"     => \$sample,
  "upstream=i"   => \$upstream,
  "downstream=i" => \$downstream
);

if ( not defined $cutoff or not defined $sample or not defined $ARGV[0] ) {
  print STDERR "[ERROR] --cutoff=90: percent of miniclusters required to have a PAS\n";
  print STDERR
    "[ERROR] --sample=id: name of the sample for which the mini-clusters are created\n";
  print STDERR "[ERROR] Usage: perl $0 --cutoff=90 --sample=id --upstream=25 --downstream=25 3pSites.PAS.tsv.gz\n";
}

if ( not -e $ARGV[0] or $ARGV[0] !~ /\.gz$/ ) {
  print STDERR "[ERROR] Could not find $ARGV[0] or file is not gzipped\n";
  exit;
}

my %data = ();
# save for each pA the line of occurrence
my $line = 0;
my %pA = ();

my $column;

open( F, "gzip -dc $ARGV[0] |" );
while (<F>) {
  chomp;
  if ( $_ =~ /^#/ ) {
    # use the header lines to infer the column for the current sample
    my @tmp_ar = split(";");
    (my $curr_sample) = $tmp_ar[1] =~ /\/([^\/]+)\.3pSites(\.highconfidence)?\.ip\.bed\.gz/;
    if( $curr_sample =~ $sample ) {
      (my $curr_col = $tmp_ar[0]) =~ s/#//;
      if( defined $column ) {
	print STDERR "[ERROR] Sample name $sample can not be associated to one input column unambiguously\n";
	exit(2);
      } else {
	$column = $curr_col;
      }
    }
    next;
  }

  # consistency check: $column has to be defined at this point
  if (not defined $column) {
    print STDERR "[ERROR] Sample $sample could not be matched with any entry available in the input data set\n";
    exit(2);
  }

  $line++;
  my @F = split(/\t/);

  #define chr:strand as key to ease access while cluster building
  my $key = "$F[0]:$F[1]";
  my $id = "$key:$F[2]";
  $pA{$id} = [ $line, $F[$column] ];
  next if( $F[$column] == 0);
  $data{$key} = [] if ( not defined $data{$key} );

  # save whether this genomic position has a PAS
  my $pas;
  if ( $F[ $#F - 1 ] eq 'NA' ) {
    $pas = 0;
  } else {
    $pas = 1;
  }

  # save for each pos: chr, strand, pos, PAS (yes/no), readCount
  push @{ $data{$key} }, [ $F[0], $F[1], $F[2], $pas, $F[$column] ];

}
close(F);

# create and save the clusters with their summed up reads
my @CLUSTERS = ();

foreach my $key ( sort keys %data ) {
  my @dataL = sort { $b->[4] <=> $a->[4] } @{ $data{$key} };

  # prepare a hash for easy access
  my %sitesHash = ();
  for ( my $i = 0 ; $i <= $#dataL ; $i++ ) {

    # save the pos in the array of each genomic pos
    $sitesHash{ $dataL[$i]->[2] } = $i;
  }

  for ( my $i = 0 ; $i <= $#dataL ; $i++ ) {
    next if ( not defined $dataL[$i] );

    my @cluster = ();
    my $start   = 0;
    my $end     = 0;
    if ( $dataL[$i]->[1] eq '+' ) {
      $start = $dataL[$i]->[2] - $upstream;
      $end   = $dataL[$i]->[2] + $downstream;
    } else {
      $start = $dataL[$i]->[2] - $downstream;
      $end   = $dataL[$i]->[2] + $upstream;
    }

    foreach my $j ( $start .. $end ) {
      if ( defined $sitesHash{$j} ) {
        my $k = $sitesHash{$j};

        # dataL[$k] contains the entry for the genomic pos $j
        next if ( not defined $dataL[$k] );
        my @deepcopy = ();
        foreach my $dc ( 0 .. $#{ $dataL[$k] } ) {
          push @deepcopy, $dataL[$k]->[$dc];
        }

        # @cluster stores all sites in the region defined as belonging to that cluster
        push @cluster, \@deepcopy;
        $dataL[$k] = undef;
      }
    }

    # get min and max cooridnates and length
    @cluster = sort { $a->[2] <=> $b->[2] } @cluster;
    my $min = $cluster[0]->[2];
    my $max = $cluster[$#cluster]->[2];
    my $l   = $max - $min + 1;

    # check if cluster has an PAS annotated
    my @tmp;
    my $pas   = 0;
    my $reads = 0;
    foreach my $site (@cluster) {

      # save positions belonging to that cluster
      push @tmp, $site->[2];
      $pas = 1 if ( $site->[3] == 1 );
      $reads += $site->[4];
    }
    push @CLUSTERS, [ $key, $reads, $pas, @tmp ];
  }
}

# define one percent of all clusters
# after the top 1% of the clusters we check the 
# fraction of clusters with PAS the first time
my $nr    = scalar @CLUSTERS;
my $firstPercent = int( $nr / 100 );

my $cnt    = 0;
my $cntPAS = 0;
my $readCutoff;
my $deleteAll = 0;
my $percent;

# do a first iteration over the clusters to find the cutoff
@CLUSTERS = sort{$b->[1] <=> $a->[1]} @CLUSTERS;
foreach my $idx ( 0 .. $#CLUSTERS ) {
  $cnt++;
  $cntPAS++ if ( $CLUSTERS[$idx]->[2] == 1 );
  next if ( $idx <= $firstPercent );
  my $pc = sprintf( "%.2f", $cntPAS / $cnt * 100 );
  $percent = $pc;

  # check the first entry after the first percent of the clusters
  if ( $idx == $firstPercent + 1 and $pc < $cutoff ) {
    print STDERR
      "[INFO] Sample $sample has very low quality. Not even among the top expressed 1 % of all clusters at least $cutoff % have a poly(A) signal.\n[INFO] Don't use this sample\n";
    print STDERR "[INFO] Percentage of the first 1 % : $percent\n";
    $readCutoff = $CLUSTERS[0]->[1];
    last;
  }
  if ( $pc < $cutoff ) {
    print STDERR "[INFO] Starting from read counts of ", $CLUSTERS[$idx]->[1],
      " mini clusters for sample $sample are not trustworthy\n";
    $readCutoff = $CLUSTERS[$idx]->[1];
    last;
  }
}
if( not defined $readCutoff ) {
  $readCutoff = 0;
  print STDERR "[INFO] all miniclusters can be considered\n[INFO] Last percentage of clusters with poly(A) signal(s): $percent\n";
}

# do a second iteration to adjust the genomic positions that should be treated as zero read counts
# for info reason: count number of deleted clusters
my $deleted_clusters = 0;
foreach my $idx ( 0 .. $#CLUSTERS ) {
  next if ( $CLUSTERS[$idx]->[1] > $readCutoff );
  $deleted_clusters++;
  my $key = $CLUSTERS[$idx]->[0];
  foreach my $pos ( 3 .. $#{ $CLUSTERS[$idx] } ) {
    my $id = $key.":".$CLUSTERS[$idx]->[$pos];
    $pA{$id}->[1] = 0;
  }
}

print STDERR "[INFO] Results based on preliminary clusters: $deleted_clusters out of $nr clusters (",
  sprintf("%.2f", $deleted_clusters / $nr * 100), " %) were deleted\n";

# output all genomic positions
foreach my $site (sort{ $pA{$a}->[0] <=> $pA{$b}->[0] } keys %pA) {
  my @F = split(":", $site);
  print join("\t", @F),"\t", $pA{$site}->[1],"\n";
}

#$mu->record('after script');
##print evaluation
#$mu->dump();##
