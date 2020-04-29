use strict;
use warnings;
use Getopt::Long;

my $minDistToPAS;
my $maxClusterSize;

GetOptions(
  "minDistToPAS=i" => \$minDistToPAS,
  "maxsize=i"      => \$maxClusterSize
);

if( not defined $minDistToPAS or not defined $maxClusterSize ) {
  print STDERR "[ERROR] Usage: perl $0 --minDistToPAS=15 --maxsize=25 clusters.primary.tsv.gz\n\n";
  exit(2);
}

if( $ARGV[0] !~ /gz$/) {
  print STDERR "[ERROR] input file does not seem to be in gzipped format ($ARGV[0] should end with .gz)\n\n";
  exit(2);
}

### ATTENTION ###
# This script expects the input file to be formatted
# according to ag-generate-clusters.pl from the A-seq-processing pipeline
# -> all data is accessed with hard coded indices
#make the meaning of indices explicit
#chr    strand  start   end     rep     total_tpm       repSite_signals [array of sampleIDs]
my %decode = ('chr' => 0, 'strand' => 1, 'start' => 2, 'end' => 3, 'rep' => 4, 'total_tpm' => 5, 'repSite_signals' => 6);

my $header;
my %clusters = ();
open( F, "gzip -dc $ARGV[0] |" ) or die "Can't open pipe to $ARGV[0]\n";
while (<F>) {
  chomp;
  if ( $_ =~ /^#/ ) {
    $header = $_;
    next;
  }
  my @F = split(/\t/);
  #save the PAS clusters in a dictionary indexed by chromosome:strand
  $clusters{"$F[0]:$F[1]"} = [] if ( not defined $clusters{"$F[0]:$F[1]"} );
  push @{ $clusters{"$F[0]:$F[1]"} }, \@F;
}
close(F);

# print header line
print "$header\n";

my $cnt_merged      = 0;
my $cnt_merged2     = 0;
my $cnt_merged3     = 0;
my $cnt_merged4     = 0;
my $cnt_merged5     = 0;
my $cnt_NA_deleted  = 0;
my $cnt_IP_deleted  = 0;
my $cnt_IP_resolved = 0;
my $clusterSelfHelp = 0;

# save the length of clusters after merging of clusters with the same poly(A) signals is done
my @lengthPASmergedCL = ();

foreach my $region ( sort keys %clusters ) {
    my @data = @{ $clusters{$region} };
    #make sure that the polyA site clusters are sorted ascendingly by the position of the representative site
    @data = sort {$a->[$decode{'rep'}] <=> $b->[$decode{'rep'}]} @data;

    # note: entries in @data are ordered in increasing
    # genomic position due to the format of the input file clusters.tsv
    
    #generate a dictionary indexed by the position in the data array, with the value being the total_tpm
    #THIS LOOP WAS BUGGY, BUT THE DATA STRUCTURE IT CREATES IS NOT USED
    my %order = ();
    for ( my $i = 0 ; $i <= $#data ; $i++ ) {
	$order{$i} = $data[$i]->[$decode{'rep'}];
    }
    
    # initially, process internal priming candidates
    ( undef, my $strand ) = split( ":", $region );
    if( $strand eq "-" ) {
	# reverse data array to iterate over it 5' -> 3'
	@data = reverse(@data);
    }
    
    for ( my $i = 0 ; $i <= $#data ; $i++ ) {
	my $repSite = $data[$i]->[4];
	# get the poly(A) signal (PAS) string
	my @rep = split( /;/, $data[$i]->[$decode{'repSite_signals' }] );
	
	#skip site if it is not marked as internal priming
	next if ( $rep[$#rep] !~ /^IP/ );

	# save all PAS-positions in the current IP-candidate-cluster
	my %ipcPAS = ();
	foreach my $entry (@rep) {
	    if ( $entry =~ m/(.*@)(\d+)$/ ) {
		$ipcPAS{$2}++;
	    }
	}

	# find next downstream cluster that is not flagged as internal priming
	my $notFound = 1;
	#start from the current cluster
	my $k     = $i;
	while ($notFound) {
	    $k++;
	    if ( $k <= $#data ) {
		if ( $data[$k]->[$decode{'repSite_signals'}] !~ /IP/ ) {
		    $notFound = 0;
		}
	    } else {
		$k     = undef;
		$notFound = 0;
	    }
	}
	if ( defined $k ) {
	    
	    # consistency check
	    if( not defined $data[$k]->[$decode{'repSite_signals'}] ) {
		print STDERR "[ERROR] No valid downstream cluster inferred\n";
		print STDERR "[ERROR] assumed downstream cluster:\n", join(",", @{$data[$k] }), "\n[ERROR] k: $k\n";
		exit(2);
	    }
	    
	    # $k is set to the downstream cluster with which the current cluster is compared
	    my $noMergePossible = 0;
	    my $sharedPAS         = 0;
	    #check all the PASes of the downstream cluster
	    my @repDS           = split( /;/, $data[$k]->[$decode{'repSite_signals'}] );
	    foreach my $entry (@repDS) {
		if ( $entry eq 'NA' ) {
		    #no PAS available, no possibility to share signals, no merging possible
		    $noMergePossible = 1;
		    last;
		}
		
		if ( $entry =~ m/(.*@)(\d+)$/ ) {
		    #downstream, non-IP cluster has some bona fide PASes
		    # check whether the one at this particular entry is also a PAS for the current cluster
		    my $entryPAS = $2;
		    if ( defined $ipcPAS{$entryPAS} ) {
			$sharedPAS = 1;
		    }
		    
		    if ( ( ($strand eq "+" and $entryPAS > $repSite ) or
			   ($strand eq "-" and $entryPAS < $repSite ) ) and not defined $ipcPAS{$entryPAS} ) {
			#the PAS of the downstream cluster is downstream of the representative cleavage site of the putative IP cluster
			#so it makes no sense to merge the clusters
			$noMergePossible = 1;
			last;
		    }
		}
	    }
	
	# check whether merging is possible
	# the non-IP cluster's PASes are all upstream of the representative site of putative IP cluster
	# and at least one PAS is shared
	if ( not $noMergePossible and $sharedPAS ) {
		
	    # construct the merged cluster
	    my $new_cl = _merge_IP( $data[$i], $data[$k] );
	    $cnt_IP_resolved++;
	    $data[$k] = $new_cl;
	    $data[$i] = undef;
	    }
	}

    #WHAT HAPPENS TO POTENTIAL IP CLUSTERS THAT WERE STILL IN BETWEEN THE IP CANDIDATE AND THE FIRST NON-IP CLUSTER DOWNSTREAM?

    # next case: cluster still exists and no merge was done
    if ( defined $data[$i] ) {
	    
	    # check whether the distance of the
	    # representative site to the closest PAS
	    # is bigger than x nt (standard x: $minDistToPAS = 15)
	    
	    ### Note for Nina: here, the strand-specific treatment was missing, causing buggy output
	    my @positions = sort { $b <=> $a } keys %ipcPAS;
	    if ( $strand eq "+" && abs( $repSite - $positions[0] ) >= $minDistToPAS ) {
		
		# cluster remains in the atlas, IP candidate note is removed
		my @tmp = @rep[ 0 .. $#rep - 1 ];
		$data[$i]->[$decode{'repSite_signals'}] = join( ";", @tmp );
		
		$clusterSelfHelp++;
		
	    } elsif ( $strand eq "-" && abs( $repSite - $positions[$#positions] ) >= $minDistToPAS ) {
		
		# cluster remains in the atlas, IP candidate note is removed
		my @tmp = @rep[ 0 .. $#rep - 1 ];
		$data[$i]->[$decode{'repSite_signals'}] = join( ";", @tmp );
		
	    } else {
		
		# cluster is discarded
		$data[$i] = undef;
		$cnt_IP_deleted++;
	    }
	}
    }
    
    # if strand is "-", change back @data entries to increasing genomic positions
    if( $strand eq "-" ) {
	@data = reverse(@data);
    }
    
    #-----------------------------------------------------------------------------
    # Finished IP-processing; next: same PAS merging
    #-----------------------------------------------------------------------------
    
    # (downstream cluster needs to have the same signals
    # as the cluster upstream of it )
    
    for ( my $i = 0 ; $i <= $#data - 1 ; $i++ ) {
	next if ( not defined $data[$i] );
	my $merge_flag = 0;
	my $k;
	my $tmp    = $i + 1;
	my $search = 1;
	while ( $search and $tmp <= $#data ) {
	    if ( defined $data[$tmp] ) {
		$k      = $tmp;
		$search = 0;
	    }
	    $tmp++;
	}
	if ( defined $k ) {
	    if ( not defined $data[$k]->[$decode{'rep'}] ) {
		print STDERR "[ERROR] data for $k inconsistent; no representative site found\n[ERROR] ", join(",", @{ $data[$k] }), "\n";
		exit(2);
	    }
	    
	    if ( $data[$i]->[$decode{'repSite_signals'}] ne 'NA' and $data[$k]->[$decode{'repSite_signals'}] ne 'NA' ) {
		# check if the representative sites are driven by the same signals
		# both clusters should have annotated PASes
		# @repA is marked as upstream, @repB as downstream cluster
		my @repA = ();
		my @repB = ();
		if ( $strand eq "+" ) {
		    @repA = split( /;/, $data[$i]->[$decode{'repSite_signals'}] );
		    @repB = split( /;/, $data[$k]->[$decode{'repSite_signals'}] );
		} elsif ( $strand eq "-" ) {
		    @repA = split( /;/, $data[$k]->[$decode{'repSite_signals'}] );
		    @repB = split( /;/, $data[$i]->[$decode{'repSite_signals'}] );
		} else {
		    print STDERR "[ERROR] Could not infer strand information from $region\n";
		    exit(2);
		}
		
		# save the PAS for the upstream cluster and merge the two clusters if the one downstream
		# has only PASes that were already found for the upstream cluster
		my %upstreamPAS         = ();
		my $PAS_present = 0;
		my $flag        = 1;
		foreach my $entry (@repA) {
		    if ( $entry =~ m/(.*@)(\d+)$/ ) {
			$upstreamPAS{$2}++;
		    }
		}
		# iterate over the PASes of the downstream cluster
		foreach my $entry (@repB) {
		    $PAS_present = 1;
		    if ( $entry =~ m/(.*@)(\d+)$/ ) {
			if ( not defined $upstreamPAS{$2} ) {
			    # donwstream cluster has a PAS not present in the upstream cluster
			    $flag = 0;
			    last;
			}
		    }
		}
		#PAS_present flag is superfluous, since both clusters should have annotated PASes
		$merge_flag = 1 if ( $flag == 1 and $PAS_present == 1 );
	    }
	    
	    if ( $merge_flag == 1 ) {
		$cnt_merged++;
		my $new_cl = _merge( $data[$i], $data[$k] );
		$data[$k] = $new_cl;
		$data[$i] = undef; #upstream cluster is removed
		
		# adjust array with cluster sizes for merged clusters
		#NOTE: merging returns a cluster with coordinates sorted in increasing order, irrespective of the strand
		$lengthPASmergedCL[$i] = 0;
		$lengthPASmergedCL[$k] = ($data[$k]->[$decode{'end'}] - $data[$k]->[$decode{'start'}]) + 1;
	    }
	}
    }
    
    #-----------------------------------------------------------------------------
    # Finished same PAS merging; next: merge if two clusters are below maxsize
    #-----------------------------------------------------------------------------
    
    # next check if the max cluster length is not yet reached
    for ( my $i = 0 ; $i <= $#data - 1 ; $i++ ) {
	next if ( not defined $data[$i] );
	my $nexti;
	my $tmp    = $i + 1;
	my $search = 1;
	while ( $search and $tmp <= $#data ) {
	    if ( defined $data[$tmp] ) {
		$nexti  = $tmp;
		$search = 0;
	    }
	    $tmp++;
	}
	if ( defined $nexti ) {
	    
	    # check if the resulting cluster is not growing too big
	    my $newL       = $data[$nexti]->[$decode{'end'}] - $data[$i]->[$decode{'start'}] + 1;
	    my $merge_flag = 0;
	    
	    # maximum length of cluster
	    if ( $newL <= $maxClusterSize ) {
		$merge_flag = 1;
		
		### Note for Nina: The following part is debatable: poly(A) clusters are
		# merged when they share the same poly(A) signals or when their maximum size is
		# below a given threshold. Depending on whether the next session is outcommented
		# or not, poly(A) clusters with different poly(A) signals are still merged when their
		# combined length is below the maximum
		
		# it is done now is it was decided by myself 2015_10_24
		# comment from this date:
		## since we use a maximum cluster length in this step
		## it seems to be more consistent to also merge
		## clusters that have independent PAS but have a 
		## combined cluster length of <= $maxClusterSize
		## -> therefore, I outcommented the region
		
		#--> from here
		
		# # compare to PASs of both clusters
		# # prevent merging, if the downstream cluster has its own PAS
		# my @repUp = ();
		# my @repDown = ();
		# my %PAStmp = ();
		# if($strand eq "+") {
		#   @repUp = split( /;/, $data[$i]->[6] );
		#   @repDown = split( /;/, $data[$nexti]->[6] );
		# } elsif( $strand eq "-") {
		#   @repUp = split( /;/, $data[$nexti]->[6] );
		#   @repDown = split( /;/, $data[$i]->[6] );
		# } else {
		#   print STDERR "[ERROR] Could not infer strand information from $region\n";
		#   exit;
		# }
		# foreach my $entry ( @repUp ) {
		#   if ( $entry =~ m/(.*@)(\d+)$/ ) {
		#     $PAStmp{$2}++;
		#   }
		# }
		# foreach my $entry ( @repDown ) {
		#   if ( $entry =~ m/(.*@)(\d+)$/ ) {
		#     if (not defined $PAStmp{ $2 } and $entry ne "NA") {
		#       # the downstream cluster has a PAS
		#       # not defined for the upstream cluster
		#       $merge_flag = 0;
		#     }
		#   }
		# }
	    }
	    
	    # <-- till here
	    
	    if ( $merge_flag == 1 ) {
		$cnt_merged2++;
		my $new_cl = _merge( $data[$i], $data[$nexti] );
		$data[$nexti] = $new_cl;
		$data[$i]     = undef;
	    }
	}
    }
    
    # Note for Nina: the following section was outcommented before
    # in the new version of the script and is active now again
    
    # process clusters without annotated PASes (NA clusters) that are close together
    my @naCluster      = ();
    my @singleClusters = ();
    for ( my $i = 0 ; $i <= $#data - 1 ; $i++ ) {
	next if ( not defined $data[$i] );
	
	if ( $data[$i]->[$decode{'repSite_signals' }] ne "NA" ) {
	    
	    # is there a cluster with no PAS to merge this cluster to?
	    if ( $#naCluster > -1 ) {
		if ( $naCluster[$decode{'end'}] - $naCluster[$decode{'start'}] <= ( $maxClusterSize * 2 ) ) {
		    
		    # merge all consecutive clusters without annotated PAS
		    if ( $#singleClusters > 0 ) {
			foreach my $entry ( 0 .. $#singleClusters - 1 ) {
			    my $newone =
				_merge( $data[ $singleClusters[$entry] ], $data[ $singleClusters[ $entry + 1 ] ] );
			    $data[ $singleClusters[$entry] ] = undef;
			    $data[ $singleClusters[ $entry + 1 ] ] = $newone;
			}
			$cnt_merged3++;
			
		    } else {
			
			# no merge necessary since there was only one previous cluster without annotated PAS
		    }
		    
		} else {
		    # distance between clusters without annotated PAS is too large
		    # discard all clusters
		    foreach my $entry ( 0 .. $#singleClusters ) {
			$data[ $singleClusters[$entry] ] = undef;
			$cnt_NA_deleted++;
		    }
		}
		@naCluster      = ();
		@singleClusters = ();
	    }
	    next;
	}
	
	# this cluster has no annotated PAS
	if ( $#naCluster == -1 ) {
	    #there is cluster without annotated PAS on the stack of consecutive such clusters
	    #put this cluster on the stack
	    @naCluster = @{ $data[$i] };
	    push @singleClusters, $i;
	    next;
	}
	
	# if we get here, then there is at least another NA cluster on the stack
	if ( $data[$i]->[$decode{'start'}] - $naCluster[$decode{'end'}] <= int($maxClusterSize/2) ) {
	    # this clusters is close to the one already on the stack; create preliminary association
	    $naCluster[$decode{'end'}] = $data[$i]->[$decode{'end'}];
	    push @singleClusters, $i;
	    
	} else {
	    # this cluster is already too far
	    # process the clusters already on the stack (which are within single linkage distance of each other)
	    if ( $naCluster[$decode{'end'}] - $naCluster[$decode{'start'}] <= ( $maxClusterSize * 2 ) ) {
		# total size of the clusters on the stack is within bounds
		# merge them
		if ( $#singleClusters > 0 ) {
		    foreach my $entry ( 0 .. $#singleClusters - 1 ) {
			my $newcluster =
			    _merge( $data[ $singleClusters[$entry] ], $data[ $singleClusters[ $entry + 1 ] ] );
			$data[ $singleClusters[$entry] ] = undef;
			$data[ $singleClusters[ $entry + 1 ] ] = $newcluster;
		    }
		    $cnt_merged3++;
		} else {
		    
		    # no merge necessary since there was a single NA cluster on the stack
		}
		
	    } else {
		# the total size of a merged cluster of all clusters currently on the stack would be too large
		# we are dealing with a broad cluster of sites without PASes, most likely spurious sites
		# discard all clusters
		# COULD BE ARGUED THAT WE SHOULD ALSO DISCARD CLUSTERS WITHOUT PAS THAT ARE NOT BROAD
		foreach my $entry ( 0 .. $#singleClusters ) {
		    $data[ $singleClusters[$entry] ] = undef;
		    $cnt_NA_deleted++;
		}
	    }
	    # we have cleared the stack, so reinitialize it with the new cluster
	    @naCluster      = @{ $data[$i] };
	    @singleClusters = ($i);
	}
    }
    
    # treat last NA cluster, if available
    if ( $#naCluster > -1 ) {
	if ( $naCluster[$decode{'end'}] - $naCluster[$decode{'start'}] <= ( $maxClusterSize * 2 ) ) {
	    
	    # merge all clusters
	    if ( $#singleClusters > 0 ) {
		foreach my $entry ( 0 .. $#singleClusters - 1 ) {
		    my $newone =
			_merge( $data[ $singleClusters[$entry] ], $data[ $singleClusters[ $entry + 1 ] ] );
		    $data[ $singleClusters[$entry] ] = undef;
		    $data[ $singleClusters[ $entry + 1 ] ] = $newone;
		}
		$cnt_merged3++;
	    } else {
		
		# no merge necessary since there was only one cluster without annotated PASes on the stack
	    }
	    
	} else {
	    
	    # the total size of a merged cluster of all clusters currently on the stack would be too large
	    # we are dealing with a broad cluster of sites without PASes, most likely spurious sites
	    # discard all clusters
	    foreach my $entry ( 0 .. $#singleClusters - 1 ) {
		$data[ $singleClusters[$entry] ] = undef;
		$cnt_NA_deleted++;
	    }
	}
	@naCluster      = ();
	@singleClusters = ();
    }

# all done, print out the final clusters   
    foreach my $entry (@data) {
	if ( defined $entry ) {
	    print join( "\t", @{$entry} ), "\n";
	}
    }
}

print STDERR
  "[INFO] $cnt_merged (same PAS) + $cnt_merged2 (below maximum length of $maxClusterSize) + $cnt_merged3 (noPAS CLs closer than 12nt to each other) entries merged.\n";
print STDERR "[INFO] $cnt_NA_deleted clusters without PAS were discarded\n";
print STDERR
  "[INFO] $cnt_IP_resolved putative internal priming clusters resolved (merged or internally solved)\n";
print STDERR "[INFO] $cnt_IP_deleted putative internal priming clusters deleted\n";
print STDERR
  "[INFO] $clusterSelfHelp IP candidates rescued since most downstream site with more than $minDistToPAS nt distance to next PAS\n";

my $CLlength_PASmergedCL = 0;
my $number               = 0;
foreach my $i ( 0 .. $#lengthPASmergedCL ) {
  if ( defined $lengthPASmergedCL[$i] and $lengthPASmergedCL[$i] > 0 ) {
    $CLlength_PASmergedCL += $lengthPASmergedCL[$i];
    $number++;
  }
}
print STDERR "[INFO] Number of clusters that emerged from same-PAS-clustering: $number\n";
print STDERR "[INFO] Average cluster length of these clusters: ",
  sprintf( "%.4f", $CLlength_PASmergedCL / $number ), "\n";

#merge two non-IP clusters
sub _merge {
    my $A = $_[0];
    my $B = $_[1];

    # determine the leader ($A entry has higher total TPM value)
    if ( $A->[$decode{'total_tpm'}] < $B->[$decode{'total_tpm'}] ) {
	$A = $_[1];
	$B = $_[0];
    }
    
    # get chr and strand
    my @C = ( $A->[$decode{'chr'}], $A->[$decode{'strand'}] );
    
    # define min start and max end
    my $minS = ( $A->[$decode{'start'}] < $B->[$decode{'start'}] ) ? $A->[$decode{'start'}] : $B->[$decode{'start'}];
    my $maxE = ( $A->[$decode{'end'}] > $B->[$decode{'end'}] ) ? $A->[$decode{'end'}] : $B->[$decode{'end'}];
    push @C, $minS;
    push @C, $maxE;
    
    # take representative site from leader cluster
    push @C, $A->[$decode{'rep'}];
    
    # sum up tpm values
    push @C, $A->[$decode{'total_tpm'}] + $B->[$decode{'total_tpm'}];
    
    # take the PAS string from $A
    push @C, $A->[$decode{'repSite_signals' }];
    
    # for each sample: sum up the expression of $A and $B
    #per-sample expression values are stored after the polyA signals
    foreach my $i ( $decode{'repSite_signals' }+1 .. $#{$A} ) {
	push @C, $A->[$i] + $B->[$i];
    }
    
    return \@C;
}

#merge an IP cluster into a non-IP cluster
sub _merge_IP {

    # $A is set as leader cluster
    my $A = $_[1];
    my $B = $_[0];

    # get chr and strand
    my @C = ( $A->[$decode{'chr'}], $A->[$decode{'strand'}] );
 
    # define min start and max end
    my $minS = ( $A->[$decode{'start'}] < $B->[$decode{'start'}] ) ? $A->[$decode{'start'}] : $B->[$decode{'start'}];
    my $maxE = ( $A->[$decode{'end'}] > $B->[$decode{'end'}] ) ? $A->[$decode{'end'}] : $B->[$decode{'end'}];
    push @C, $minS;
    push @C, $maxE;

    # take representative site from leader cluster
    push @C, $A->[$decode{'rep'}];

    # sum up tpm values
    push @C, $A->[$decode{'total_tpm'}] + $B->[$decode{'total_tpm'}];

    # take the PAS string from the leader cluster
    push @C, $A->[$decode{'repSite_signals' }];

    # for each sample: sum up the expression of $A and $B
    #per-sample expression values are stored after the polyA signals
    foreach my $i ( $decode{'repSite_signals' }+1 .. $#{$A} ) {
	push @C, $A->[$i] + $B->[$i];
    }
    
    return \@C;
}
