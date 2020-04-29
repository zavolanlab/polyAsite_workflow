use strict;
use warnings;
use Getopt::Long;

# Define/initialize user-editable variables
my $upstream;                   # From CLI: number of nucleotides to be considered upstream of high read support sites when determining clusters; maximum cluster size is given by `$upstream` + `$downstream` + 1
my $downstream;                 # From CLI: number of nucleotides to be considered downstream of high read support sites when determining clusters; maximum cluster size is given by `$upstream` + `$downstream` + 1
my $first_cnt_col       = 4;    # 1-based index of first count column (e.g., if counts start in fourth column, put 4)
my $minTPM              = 0;    # Minimum total expression of a cluster for it to be reported

# Define/initialize other variables
my $n_clusters          = 0;    # Number of valid clusters
my $ip_candidates       = 0;    # Counter for clusters that are candidates for internal priming
my $total_TPM           = 0;    # Total expression of all clusters
my $TPM_clusters_gt_1nt = 0;    # Total expression of all clusters extending over more than 1 nucleotide
my %data                = ();   # Hash of arrays of arrays, with "chromosome:strand" as keys; arrays hold counts/expression values etc., as well as genomic info required for clustering per row of data (i.e., genomic position)
my @LIBSIZE             = ();   # Array of library sizes in order of 0-based column index in input file header, starting with index 3 (first "counts" column)
my @SMPL_ID             = ();   # Array of sample identifiers in order of 0-based column index in input file header, starting with index 0
my $N                   = 0;    # Current data row of input file (not counting header lines)
my $nr_of_samples;              # Number of samples (and count data columns) in input file; corresponds to length of non-null values in @LIBSIZE and @SMPL_ID

# Get command line options
GetOptions(
    "upstream=i"   => \$upstream,
    "downstream=i" => \$downstream
    );

# Valid inputs & CLI options
if ( not defined $ARGV[0] ) {
    print STDERR "[ERROR] No input file.\n";
    exit(2);
}
if ( not defined $upstream or not defined $downstream ) {
    print STDERR "[ERROR] Usage: perl $0 --upstream=12 --downstream=12 3pSites.PAS.noBG.tsv.gz\n";
    exit(2);
}

# Open and read entire input file (header and body), set `$nr_of_samples`,
# `@SMPL_ID`, `@LIBSIZE` and add necessary data to hash `%data`
# Note: Read counts for each sample are listed in the columns starting with
# `$first_cnt_col`
open( F, "gunzip -c $ARGV[0] |" );
while (<F>) {

    # Define/initialize local variables
    my $OVERALL_TPM = 0;   # Mean TPM across all samples
    my $n_supported = 0;   # Number of samples with TPM > 0
    my @tmp         = ();  # Array holding normalized, formatted counts per sample (TPM; "tags per million")

    # Process line
    chomp;
    # Get sample ID and library size from header line
    # example format: #3;/scicore/home/zavolan/herrmchr/PolyASite/Data/samples_20190730_cl_avg_tpm/GSM1614163/GRCh38-96/GSM1614163.3pSites.ip.bed.gz;605571
    # as long as there are header lines (starting with comment) the block up to next is executed
    if ( $_ =~ m/^#(.+)/ ) {
	my @F = split( /;/, $1 ); #example: (3, {path}, 605571)
	my @G = split( /\//, $F[1]); #parse the path
	my @H = split( /\./, $G[-1]); #parse the file name to remove extensions
	push @SMPL_ID, $H[0]; #keep track only of sample ID
	push @LIBSIZE, $F[2]; #last number in header line is sample size
	next;
    }
    # We only get here when the header lines exhausted
    # Split body line by whitespace
    my @F = split(/[\t\s]+/);
    #keep track of the number of polyA sites
    $N++;

    # Get numer of samples
    if ( not defined $nr_of_samples ) { $nr_of_samples = scalar @SMPL_ID };

    # Process the read counts per site from each sample
    foreach my $idx ( 0 .. $nr_of_samples - 1 ) {
	#extract count
	my $cnt = $F[$idx + $first_cnt_col - 1];
	#normalize to TPM
	my $tpm = $cnt / $LIBSIZE[$idx] * 1000000;
	#save the formatted value of the normalized TPM
	push @tmp, sprintf( "%.6f", $tpm);
	#keep track of the number of samples where the site was observed
	$n_supported++ if ( $tpm > 0 );
	#keep track of cumulative TPM to be able to calculate the average
	$OVERALL_TPM += $tpm;
    }
    # Take mean of total TPM per row of data
    $OVERALL_TPM = $OVERALL_TPM / $nr_of_samples;

    # Create hash "%data" with "chromosome:strand" as keys; initiate with empty list
    my $key = "$F[0]:$F[1]";
    $data{$key} = [] if ( not defined $data{$key} );

    # Save chr, strand, pos, overall tpm, number of supporting samples,
    # the poly(A) signals and the tpm values per sample
    push @{ $data{$key} }, [ $F[0], $F[1], $F[2], $OVERALL_TPM, $n_supported, $F[ $#F - 1 ], @tmp ];
}
close(F);

# Print log info & header line
print STDERR "[SAMPLES] @SMPL_ID\n";
print STDERR "[INFO] $N sites read.\n";
print "#chr\tstrand\tstart\tend\trep\ttotal_tpm\trepSite_signals\t", join("\t", @SMPL_ID), "\n";

#data hash is indexed by "chr:strand"
#values are: chr, strand, pos, overall tpm, number of supporting samples, poly(A) signals, tpm values per sample
#we'll save these fields in a hash
my %decode = ('chr' => 0, 'strand' => 1, 'position' => 2, 'overallTPM' => 3, 'supportingSamples' => 4, 'polyAsignals' => 5);

# Sort sites per "chromsome:strand"
foreach my $key ( sort keys %data ) {

  # Print log message
  print STDERR "[INFO] Processing '$key'.\n";

  # Define/initialize local variables
  my @OUTPUT = ();  # Array to hold values to be written to STDOUT

  # Sort sites by `$n_supported`, then by `$OVERALL_TPM`
  my @dataL = sort { $b->[$decode{'supportingSamples'}] <=> $a->[$decode{'supportingSamples'}] || $b->[$decode{'overallTPM'}] <=> $a->[$decode{'overallTPM'}] } @{ $data{$key} };

  # For easy access, prepare hash with genomic positions as keys & `@dataL`
  # indices as values
  my %sitesHash = ();
  for ( my $i = 0 ; $i <= $#dataL ; $i++ ) {
    $sitesHash{ $dataL[$i]->[$decode{'position'}] } = $i;
  }

  # Cluster sites in region defined by `$upstream` and `$downstream` around
  # position of current site, starting with the most supported/"expressed"
  # sites
  for ( my $i = 0 ; $i <= $#dataL ; $i++ ) {

    # Define/initialize local variables
    my @cluster     = ();  # Array of arrays of clustered entries in (a sorted copy of) `%data`
    my $start       = 0;   # Start position (lowest number) of region/window around current site
    my $end         = 0;   # End position (higest number) of region/window around current site
    my %rep         = ();  # Hash of positions (keys) and times (values; integers) that that position had the highest normalized read count across all samples
    my $N_supported = 0;   # Number of samples with a non-zero read count in a given cluster
    my @sumTPM      = ();  # Array of sums of normalized read counts of each sample in a given cluster

    # Skip if data is not available
    next if ( not defined $dataL[$i] );  # This shouldn't happen as an error would be thrown in the loop above

    # Define region/window around position of current site
    if ( $dataL[$i]->[1] eq "+" ) {
      $start = $dataL[$i]->[$decode{'position'}] - $upstream;
      $end   = $dataL[$i]->[$decode{'position'}] + $downstream;
    } else {
      $start = $dataL[$i]->[$decode{'position'}] - $downstream;
      $end   = $dataL[$i]->[$decode{'position'}] + $upstream;
    }

    # Iterate over each position of region/window around current site, and if 
    # a site has been defined for that position, add its data to `@cluster`
    foreach my $j ( $start .. $end ) {

      # Check if site has been defined
      if ( defined $sitesHash{$j} ) {

        # Get `@dataL` index of position
        my $k = $sitesHash{$j};

        # Skip if data for position is not available
        next if ( not defined $dataL[$k] );

        # Add data to cluster
        my @deepcopy = ();
        foreach my $dc ( 0 .. $#{ $dataL[$k] } ) {
          push @deepcopy, $dataL[$k]->[$dc];
        }
        push @cluster, \@deepcopy;

        # To avoid double counting, delete data that has been added to cluster
        $dataL[$k] = undef;
      }

    }

    # Get strand, min/max positions and length of cluster
    @cluster = sort { $a->[$decode{'position'}] <=> $b->[$decode{'position'}] } @cluster;
    my $strand = $cluster[0]->[1];
    my $min = $cluster[0]->[$decode{'position'}];
    my $max = $cluster[$#cluster]->[$decode{'position'}];
    my $len = $max - $min + 1;

    ## Get representative site
    # Define representative site as the site that most often had the most
    # read support across all samples; in case of equality, the most
    # downstream site is chosen

    # Iterate over samples
    foreach my $idx ( 6 .. $#{ $cluster[0] } ) {

      # Sort sites in cluster by normalized read count in current sample
      my @tmp = sort { $b->[$idx] <=> $a->[$idx] } @cluster;

      # If the current sample has at least one non-zero read count for any site
      # in the current cluster, increase by 1 the counter `$N_supported` as well
      # as the counter `$rep{$pos}` for all positions `$pos` that have the
      # highest read count
      my $max_count = undef;
      my $sum_tmp = 0;
      foreach my $curr_site (@tmp) {

	  # Break out of loop if site has no read support for current sample,
	  # else increase counter `$N_supported`
	  # Sites are sorted in decreasing order of counts, so once we hit 0 counts, we can get out of the loop
        if ( $curr_site->[$idx] > 0 ) { $N_supported++; } else { last; }

        # Add normalized read count of current site to sum of normalized read
        # count per sample
        $sum_tmp += $curr_site->[$idx];

        # Set highest read count (determined by first site after sorting!)
	# This works because the sites are sorted in decreasing order of counts
        if ( not defined $max_count ) { $max_count = $curr_site->[$idx]; }

        # Increase `$rep{$pos}` counter for current position if current site
        # has highest read count
        if ( $curr_site->[$idx] == $max_count ) { $rep{ $curr_site->[$decode{'position'}] }++; }

      }
      # Add read count sum to array `@sumTPM`
      push @sumTPM, sprintf( "%.6f", $sum_tmp);

    }

    # Sort keys (i.e., positions) of counter `%rep` by read support first,
    # then by position; position sorting is dependent on strand, from the
    # most downstream to the most upstream position
    my @R;
    if ( $strand eq "+" ) {
      @R = sort { $rep{$b} <=> $rep{$a} || $b <=> $a } keys %rep;
    } else {
      @R = sort { $rep{$b} <=> $rep{$a} || $a <=> $b } keys %rep;
    }
    # Representative site is now defined by the first element of the sorted
    # array
    my $rep_pos = $R[0];

    # Calculate the fraction of samples that support the representative site
    # and the cluster
    my $rep_support = sprintf( "%.2f", $rep{ $rep_pos } / $nr_of_samples );
    $N_supported = sprintf( "%.2f", $N_supported / $nr_of_samples );

    # Find internal priming (IP) candidates within cluster and get poly(A)
    # signal (PAS) string for representative site
    # PAS string example:
    # AATAAA@-29@125193;AATATA@-24@125198
    # i.e. poly(A) signal, position of PAS start relative to cleavage site, absolute position in the genome 
    # Here, IP candidates are defined as 3p/cleavage sites overlapping directly
    # with the PAS, which, due to its A/T richness, may have prompted the priming
    my $rep_PAS = 'NA';
    my $ipCandidate = 0;
    foreach my $site ( 0 .. $#cluster ) {

      # Extract PAS string
      my $pasString = $cluster[$site]->[$decode{'polyAsignals'}];

      # Skip entry if no PAS is available
      next if ( $pasString eq 'NA' );

      # Get poly(A) signal (PAS) string for representative site, if available
      if ( $cluster[$site]->[2] == $rep_pos ) {
        $rep_PAS = $cluster[$site]->[$decode{'polyAsignals'}];
      }

      # Do only if not already marked as IP candidate
      if ( ! $ipCandidate ) {

        # Split multiple PAS
        my @signals = split( ";", $pasString );

        # Iterate over PAS
        foreach my $sig (@signals) {

          # Extract relative position
          ( my $pas, my $pos ) = split( "@", $sig );

          # Mark cluster as IP candidate if relative position is within 6 nts
          # of PAS
          if ( $pos >=  -length($pas)) {
            $ipCandidate = 1;
            last;
          }

        }

      }
    
    }

    # Update representative sites & PAS string for IP candidates
    # Note: For IP candidates, the most downstream site is defined as the
    # representative site
    if ($ipCandidate) {

      # Increase IP candidate counter
      $ip_candidates++;

      # Sort cluster by genome position such that most downstream site
      # (strand-dependent!) is the first element in array '@cluster'
      if ( $strand eq "+" ) {
        @cluster = sort { $b->[2] <=> $a->[2] } @cluster
      } else {
        @cluster = sort { $a->[2] <=> $b->[2] } @cluster
      }

      # Update representative site and related info
      $rep_pos     = $cluster[0]->[2];
      $rep_PAS     = $cluster[0]->[$decode{'polyAsignals'}];
      $rep_support = 0;
      $N_supported = 0;

      # Add "IPCandidate" to the PAS string
      $rep_PAS .= ";IPCandidate";

    }

    # Calculate total expression for current cluster
    my $cluster_tpm = 0;
    foreach my $site ( @cluster ) {
      $cluster_tpm += $site->[$decode{'overallTPM'}];
    }

    # Prepare output & stats if cluster expression is above specified threshold
    if ( $cluster_tpm >= $minTPM ) {

      # Increase counter of valid clusters
      $n_clusters++;

      # Add total expression for current cluster to total expression
      # `$total_TPM` of all clusters
      $total_TPM += $cluster_tpm;

      # Add total expression for current cluster to total expression
      # `$TPM_clusters_gt_1nt` of all clusters extending over more than 1
      # nucleotide
      $TPM_clusters_gt_1nt += $cluster_tpm if ( $len > 1 );

      # Prepare output array
      push @OUTPUT, [
        $cluster[0]->[0],
        $strand, $min, $max, $rep_pos, $cluster_tpm, $rep_PAS, @sumTPM
        ];

    }

  }
  #make the meaning of indices obvious again
  my $rep_position_output = 4;
  # Print output, sorted by representative position
  @OUTPUT = sort { $a->[$rep_position_output] <=> $b->[$rep_position_output] } @OUTPUT;
  foreach my $entry (@OUTPUT) {
    print join( "\t", @{$entry} ), "\n";
  }

}

# Print stats to log
my $pc = sprintf( "%.2f", $TPM_clusters_gt_1nt / $total_TPM * 100 );
$total_TPM = sprintf( "%.2f", $total_TPM);
print STDERR
  "[INFO] Input file: $ARGV[0]\n",
  "[INFO] Clusters: $n_clusters\n",
  "[INFO] Total expression of all clusters: $total_TPM TPM\n",
  "[INFO] Percent of expression associated with clusters > 1nt: $pc%\n",
  "[INFO] Internal priming candidate clusters: $ip_candidates\n",
  "[INFO] Clusters were inferred with distances of $upstream upstream and ",
  "$downstream downstream nucleotides\n[INFO] from sites with highest ",
  "sample/read support. Clusters were required to have a\n[INFO] total ",
  "expression of at least $minTPM TPM in order to be considered.\n";
