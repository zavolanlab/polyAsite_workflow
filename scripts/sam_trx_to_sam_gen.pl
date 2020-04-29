#!/usr/bin/env perl
use lib "/scicore/home/zavolan/clipz/newClipz6/lib/perl";
#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: sam_trx_to_sam_gen.pl
### Created: Oct 4, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Requirements: Getopt::Long
#==================#
#    HEADER END    #
#==================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#

#---> PRAGMAS / PACKAGES <---#
use strict;
use warnings;
use Getopt::Long;

#---> USAGE <---#
my $usage_info = &usage;

#---> OPTIONS / ARGUMENTS <---#
my $usage = '';
my $quiet = '';
my $in_sam = '';
my $in_bed = '';
my $out_sam = '';
my $head = '';
my $no_strand_info = 0;
my $min_overlap = 0;
my $monoexonic = 0;
my $tag = '';
my $report = 0;
my $options_result = GetOptions (
	'usage|help' => \$usage,
	'quiet' => \$quiet,
	'head' => \$head,
	'no-strand-info' => \$no_strand_info,
	'min-overlap=i' => \$min_overlap,
	'include-monoexonic' => \$monoexonic,
	'tag=s' => \$tag,
	'print-report' => \$report,
	#-----------------------#
	'in=s' => \$in_sam,
	'exons=s' => \$in_bed,	
	'out=s' => \$out_sam
);

# Die if option parsing was not successful or --usage / --help was requested
die $usage_info if $usage || !$options_result;

# Die if required options are not provided
die $usage_info if !$in_sam || !$in_bed || !$out_sam; 

# Construct optional tag
$tag = "\tZZ:Z:" . $tag if $tag; 

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#

#---> STATUS MESSAGE <---#
print "Starting '$0'...\n" unless $quiet;

#---> MAIN VARIABLES <---#
my $exon_hoaoa_ref;

#---> BODY <---#

	#---> Construct hash of arrays of arrays of exons <---#
	$exon_hoaoa_ref = &exons_bed_to_hoaoa($in_bed);
#	$exon_hoaoa_ref = &exons_hoaoa_remove_single_exons($exon_hoaoa_ref) unless $monoexonic;
	$exon_hoaoa_ref = &exons_hoaoa_reorder_exons_on_minus_strand($exon_hoaoa_ref);
	$exon_hoaoa_ref = &exons_hoaoa_add_cumulative_length($exon_hoaoa_ref);
	$exon_hoaoa_ref = &exons_hoaoa_add_intron_length($exon_hoaoa_ref);

	#---> Map reads <---#
	&trx_sam_to_gen_sam($in_sam, $exon_hoaoa_ref, $out_sam);

#---> STATUS MESSAGE <---#
print "Done.\n" unless $quiet;

#---> PROGRAM EXIT <---#
exit 0;

#================#
#    MAIN END    #
#================#


#=======================#
#   SUBROUTINES START   #
#=======================#
sub usage {
### Function: Returns usage information for current script
### Accepts: n/a
### Returns: String with usage information
### Type: Specialized
'Usage: perl ./sam_trx_to_sam_gen.pl [OPTIONS] --in [FILE|SAM] --exons [FILE|BED] --out [FILE|SAM]

Description: Re-maps a SAM file resulting from aligning a library of sequencing reads against a transcriptome to genomic coordinates. All reads that do not cross an exon-exon boundary by a specified minimum overlap are discarded by default.

==================================================
Required arguments:
--in [FILE|SAM]	Input SAM file (transcript coordinates)
--exons [FILE|BED]	BED file of exons (genomic coordinates)
--out [FILE|SAM]	Output SAM file (genomic coordinates)
==================================================
Optional arguments:
--min-overlap [INT]	Minimum required overlap between read and any of the exons
--no-strand-info	Used library preparation protocol does not preserve strand information (if unset, all reads mapping to the opposite strand of annotated transcripts - i.e. the appropriate SAM flag is set - are discarded)
--include-monoexonic	Do not discard alignments against single exons (default: skip)
--head	Print SAM header
--tag [STRING]	Tag of the form "ZZ:Z:STRING" that is appended to the end of each line in the output file
--print-report	Print statistics on the number of processed, printed and discarded alignments to STDOUT
--usage|help	Show this information and die
--quiet	Shut up!

Comments:
- The script was written for SAM files produced by a recent (>= 0.1.4) version of segemehl. SAM files of other origin may or may not work. 
- Only the 0x10 bit (which informs whether a read sequence was reverse complemented for the alignment) of the FLAG field is considered; all other bits are ignored.
- CAUTION: Only marginal validation of the input file type/format performed!

Version 1.6 (2014-08-04)
Written by Alexander Kanitz on 2013-10-04
';
}
#-----------------------#
sub exons_bed_to_hoaoa {
### Function: Reads a BED file of exons and loads them into a hash (key: transcript ID) of arrays (chromosome, strand, exons) of arrays (exon genomic start, exon genomic stop)
### Accepts: 1. sorted (name, start, end) BED6 file of exons (coordinates relative to genome; 0-based, open-ended)
### Returns: Reference to hash of arrays of arrays
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $bed = shift;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Loading exons..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my %hoaoa;
	my @outer_array;
	
	#---> BODY <---#

		#---> Open input filehandle <---#
		open IN, "<", $bed or die "[ERROR] Could not open file '$bed'!\nExecution aborted\n";
		
		#---> Traverse input file line by line <---#
		while (my $line = <IN>) {		
			
			#---> Field values to array <---#
			chomp $line;
			my ($chr, $start, $end, $name, $score, $strand) = split "\t", $line;
			die "[ERROR] Input file does not look like a valid BED6 file!\nExecution aborted\n" unless $strand;	# Assert presence of strand information
			warn "[WARNING] Incompatible strand information in line $.:\n$line\nOnly '+' and '-' are allowed. Entry skipped.\n" and next unless ($strand eq "+") || ($strand eq "-");
			
			#---> Generate outer array (chromosome, strand) if no record yet exists for this transcript <---#			
			if (! exists $hoaoa{$name}) {
				undef @outer_array;
				push @outer_array, $chr, $strand;
			}
	
			#---> Verify that one transcript does not contain exons annotated on different strands <---#
			die "[ERROR] Transcript $name contains exons annotated on different strands in/around line $.:\n$line\nExecution aborted\n" if $outer_array[1] ne $strand;

			#---> Generate inner array ()exon length, genomic start / end) and push to outer array <---#
			my @inner_array;
			push @inner_array, ($start + 1), $end;	# Transform to 1-based, close-ended (SAM!)

			#---> Push inner to outer array and add outer array to hash <---#
			push @outer_array, \@inner_array;
			$hoaoa{$name} = [ @outer_array ];
		}
		
		#---> Close input filehandle <---#
		close IN;			
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Exons loaded." . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \%hoaoa;
	
}
#-----------------------#
sub exons_hoaoa_reorder_exons_on_minus_strand {
### Function: Reverses the order of exons for transcripts annotated on the Crick/minus strand from exons hoaoa generated by sub "exon_bed_to_hoaoa"
### Accepts:  Reference to hash of arrays of arrays
### Returns:  Reference to hash of arrays of arrays
### Dependencies: Subroutine "exon_bed_to_hoaoa" (Alexander Kanitz)
### Type: Specialized
	#---> PASS ARGUMENTS ---#
	my $hoaoa_ref = shift;	
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Reversing exon order for transcripts on Crick strand..." . "\n" unless $quiet;

	#---> BODY <---#

		#---> Reverse order of exons for transcripts annotated on Crick strand <---#
		foreach my $outer_array_ref ( values %$hoaoa_ref ) {
			# If transcript is annotated on Crick strand
			if ( $outer_array_ref->[1] eq "-" ) {
				# Slice out chr and strand
				my @chr_str = splice @$outer_array_ref, 0, 2;
				## Swap start and end coordinates
				foreach my $exon_array_ref (@$outer_array_ref) {
					my $old_start_new_stop = $$exon_array_ref[0];
					$$exon_array_ref[0] = $$exon_array_ref[1];
					$$exon_array_ref[1] = $old_start_new_stop;
				}
				# Reverse rest (i.e. references to exons / inner arrays)
				my @reversed_inner_array_refs = reverse @$outer_array_ref;
				# Empty original array
				@$outer_array_ref = ();
				# Rebuild array from slices
				push @$outer_array_ref, @chr_str, @reversed_inner_array_refs;
			}
		}	
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Exon order reversed." . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return $hoaoa_ref;
	
}
#-----------------------#
sub exons_hoaoa_add_cumulative_length {
### Function: Add exon length and cumulative exon length to exons hoaoa generated by sub "exon_bed_to_hoaoa" (added as third and fourth elements to inner arrays)
### Accepts: Reference to hash of arrays of arrays
### Returns: Reference to hash of arrays of arrays
### Dependencies: Subroutine "exon_bed_to_hoaoa" (Alexander Kanitz)
### Type: Specialized
	#---> PASS ARGUMENTS ---#
	my $hoaoa_ref = shift;	
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Calculating cumulative exon lengths..." . "\n" unless $quiet;

	#---> BODY <---#

		#---> Add cumulative exon length to inner arrays <---#
		## Traverse through each key $outer_array of hash %$hoaoa_ref
		foreach my $outer_array_ref ( values %$hoaoa_ref ) {
			# Initialize cumulative length variable
			my $previous = 0;
			## For each inner array reference...
			for my $exon_no ( 0 .. (scalar @$outer_array_ref - 3) ) {
				# Calculate current exon length 
				$outer_array_ref->[$exon_no + 2]->[2] = abs($outer_array_ref->[$exon_no + 2]->[1] - $outer_array_ref->[$exon_no + 2]->[0]) + 1;
				# Calculate cumulative length by adding current exon length to previous cumulative length; assign to new element in inner array
				$outer_array_ref->[$exon_no + 2]->[3] = $outer_array_ref->[$exon_no + 2]->[2] + $previous;
				# Update cumulative length
				$previous = $outer_array_ref->[$exon_no + 2]->[3];
			}
		}
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Cumulative exon lengths calculated." . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return $hoaoa_ref;
	
}
#-----------------------#
sub exons_hoaoa_add_intron_length {
### Function: Add distance between exon and previous exon to exons hoaoa generated by sub "exon_bed_to_hoaoa" (added as fifth element to inner arrays)
### Accepts: Reference to hash of arrays of arrays
### Returns: Reference to hash of arrays of arrays
### Dependencies: Subroutine "exon_bed_to_hoaoa" (Alexander Kanitz)
### Type: Specialized
	#---> PASS ARGUMENTS ---#
	my $hoaoa_ref = shift;	
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Calculating intron lengths..." . "\n" unless $quiet;

	#---> BODY <---#

		#---> Add length of previous intron to inner arrays <---#
		## Traverse through each key $outer_array of hash %$hoaoa_ref
		foreach my $outer_array_ref ( values %$hoaoa_ref ) {
			# Initialize variable that holds end coordinate of previous exon
			my $end_prev;
			## For each inner array reference...
			for my $exon_no ( 0 .. (scalar @$outer_array_ref - 3) ) {
				# Calculate intron length by subtracting the end coordinate of the previous exon from the start coordinate of the current exon (adjust offset!); set 0 if there is no previous exon; assign value to new element in inner array
				$outer_array_ref->[$exon_no + 2]->[4] = (defined $end_prev) ? abs($outer_array_ref->[$exon_no + 2]->[0] - $end_prev) - 1 : 0;
				# Update/set end coordinate for next iteration
				$end_prev = $outer_array_ref->[$exon_no + 2]->[1];
			}
		}
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Intron lengths calculated." . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return $hoaoa_ref;
	
}
#-----------------------#
sub trx_sam_to_gen_sam {
### Function: Reads transcript alignments from a SAM file and maps them to genomic coordinates
### Accepts: 1. SAM input file; 2. Hash of arrays of arrays of exon genomic coordinates generated by subroutine "exons_bed_to_hoaoa" or derivatives; 3. Filename for SAM output file
### Returns: n/a
### Dependencies: Subroutine "exons_bed_to_hoaoa" (Alexander Kanitz)
### Type: Specialized
	#---> PASS ARGUMENTS ---#
	my ($in_sam, $exons_hoaoa_ref, $out_sam) = @_;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Processing SAM file (may take long)..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my $pos;
	my $processed = 0;
	my $incomplete = 0;
	my $rc = 0;
	my $not_found = 0;
	my $no_eej = 0;
	my $below_min_overlap = 0;
	my $printed = 0;
	
	#---> BODY <---#
	
		#---> Open input and output filehandles <---#
		open IN, "<", $in_sam or die "[ERROR] Could not open file '$in_sam'!\nExecution aborted\n";
		open OUT, ">", $out_sam or die "[ERROR] Could not open file '$out_sam'!\nExecution aborted\n";
		
		#---> Process header lines (assumed at the top of the file!) <---# 
		$pos = tell IN;
		while (<IN>) {
			if ( $_ =~ /^\@[A-Za-z][A-Za-z]\s/ ) {
				print OUT if $head;		# Print header if --head switch is set
				$pos = tell IN;
			}
			else {	
				last;			
			}
		}
	
		#---> Reset filehandle position  <---#
		seek IN, $pos, 0;
		$. = $. - 1;

		#---> Traverse non-header lines <---#
		while (<IN>) {

			#---> Process SAM line <---#
			chomp;
			next if $_ eq "";
			$processed++;
			my ($QNAME, $FLAG, $RNAME, $POS, $MAPQ, $CIGAR, $RNEXT, $PNEXT, $TLEN, $SEQ, $QUAL, $REST) = split "\t", $_, 12;

			#---> Check if all (i.e. last) required values are defined <---#
			unless (defined $QUAL) {
				warn "[WARNING] Line ${.} does not look like a proper SAM line. Entry skipped.\n";
				$incomplete++;
				next;
			}
			

			#---> Extract MD tag <---#			
			my $MD;
			if ( defined $REST ) {
				$REST =~ s/\t(MD:Z:\S+)//;
				$MD = $1;
			}
			
			#---> Skip reads when alignments are reverse complemented but strand information is available <---#
			if ( &rc_bit($FLAG) == 1 && ! $no_strand_info ) {
				$rc++;
				next;
			}
			
			#---> Adjust start and calculate end position of transcript alignment <---#
			my $ref_length = 0;
			$ref_length += $_ for $CIGAR =~ /(\d+)[MDN]/g;		# Only M(ism)atches, D(eletions) and Ns contribute to the aligned part of the reference
			my $END_POS = $POS + $ref_length - 1;

			#---> Obtain genomic coordinates, lengths of intersections with exons and lengths of spanned introns <---#
			my $coords_frags_array_ref = &get_coords_and_frags($exons_hoaoa_ref, $RNAME, $POS, $END_POS, $min_overlap, $monoexonic);
			if ($coords_frags_array_ref == 1) {			# Skip alignment if previous function call returned error code 1 (i.e. aligns to a transcript or transcript region not in the lookup table)
				$not_found++;
				next; 
			}				
			if ($coords_frags_array_ref == 2) {			# Skip alignment if previous function call returned error code 2 (i.e. alignment does not intersect multiple exons)
				$no_eej++;
				next;
			}
			if ($coords_frags_array_ref == 3) {			# Skip alignment if previous function call returned error code 3 (i.e. the specified minimum overlap is not met)
				$below_min_overlap++;
				next;
			}
			my ($chr, $str, $start, @frags) = @$coords_frags_array_ref;

			#---> Evaluate strand information and update FLAG, SEQ, QUAL, CIGAR & MD if necessary <---#
			if ( $str eq "-" ) {
				$FLAG = ($FLAG == 0) ? 16 : 0;		# The Watson ("+") strand is the reference strand for the genome alignments, therefore a perfect match (FLAG = 0) to a transcript on the "+" and "-" strands should be assigned FLAG values of 0 and 16, respectively. The situation is inversed if a reverse complement match occurs (only considered when --no-strand-info is set!)
				$CIGAR = &reverse_CIGAR($CIGAR);
				$SEQ = reverse &complement($SEQ);
				$QUAL = reverse $QUAL;
				$MD = &reverse_complement_MD($MD) if defined $MD;
			};

			#---> Add introns to CIGAR string <---#
			$CIGAR = &add_introns_to_CIGAR($CIGAR, @frags) unless scalar @frags == 0;

			#---> Print entry <---#
			print OUT
					$QNAME	. "\t" .
					$FLAG	. "\t" .
					$chr	. "\t" .
					$start	. "\t" .
					$MAPQ	. "\t" .
					$CIGAR	. "\t" .
					$RNEXT	. "\t" .
					$PNEXT	. "\t" .
					$TLEN	. "\t" .
					$SEQ	. "\t" .
					$QUAL	. "\t" .
					$REST;
			print OUT $tag if $tag;
			print OUT "\t" . $MD if defined $MD;
			print OUT "\n";
			$printed++;
		}
		
		#---> Print report <---#		
		my $discarded = $incomplete + $rc + $not_found + $no_eej + $below_min_overlap;
		my $stats = "======\nREPORT\n======\n";
		$stats .= sprintf "Alignments processed: %d (%.2f%%)\n", $processed, $processed / $processed * 100;
		$stats .= sprintf "  - Converted to genome space: %d (%.2f%%)\n", $printed, $printed / $processed * 100;
		$stats .= sprintf "  - Discarded: %d (%.2f%%)\n", $discarded, $discarded / $processed * 100;
		$stats .= sprintf "    - Sequence reverse complemented: %d (%.2f%%)\n", $rc, $rc / $processed * 100;
		$stats .= sprintf "    - Not covering exon-exon junctions: %d (%.2f%%)\n", $no_eej, $no_eej / $processed * 100;
		$stats .= sprintf "    - Not meeting minimum exon overlap: %d (%.2f%%)\n", $below_min_overlap, $below_min_overlap / $processed * 100;
		$stats .= sprintf "    - Features/regions not present in exon lookup table: %d (%.2f%%)\n", $not_found, $not_found / $processed * 100;
		$stats .= sprintf "    - Incomplete/non-standard records: %d (%.2f%%)\n", $incomplete, $incomplete / $processed * 100;
		
		print STDOUT $stats if $report;

		#---> Close input and output filehandles <---#
		close OUT;
		close IN;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "SAM file processed. Output written to '$out_sam'" . "\n" unless $quiet;

}
#-----------------------#
sub get_coords_and_frags {
### Function: Gets the reference sequence (i.e. chromosome), strand, and starting position of the alignment in genomic coordinates and calculates the lengths of overlaps with intersected exons as well as the lengths of spanned introns
### Accepts: 1. Reference to hash of arrays of arrays generated by subroutine "exons_bed_to_hoaoa" or derivatives; 2./3. RNAME (i.e. transcript ID) and POS (i.e. starting position) of the transcript SAM entry, respectively; 4. end position of the transcript alignment; 5. the allowed minimum overlap
### Returns: Reference to array of 1. chromosome, 2. strand, 3. starting position of the alignment in genomic coordinates, 4. overlap with first exon (integer), (for multiple fragments, alternating: A. length of spanned intron, B. overlap with next fragment), N-1. length of spanned intron, N. overlap with last fragment
### Dependencies: Subroutine "exons_bed_to_hoaoa" (Alexander Kanitz)
### Type: Specialized
	#---> PASS ARGUMENTS ---#
	my ($exons_hoaoa_ref, $trx_id, $start_trx, $end_trx, $min_overlap, $monoexonic) = @_;
	
	#---> SUBROUTINE VARIABLES <---#
	my @return;
	my @frags;
	my $start_gen;
	my $next_index;
	my $prev_cum_length = 0;
	my $skip = 1;

	#---> BODY <---#

		#---> Get array containing properties of transcript <---#
		# Throw warning and return an error code if transcript is missing from lookup table
		warn "[WARNING] Transcript '$trx_id' not found in lookup table! Check integrity and compatibility of exon/transcript annotations used here and during mapping. Entry skipped.\n" and return 1 unless exists $exons_hoaoa_ref->{$trx_id};
        # Load entry into dedicated array
        my @trx_entry = @{$exons_hoaoa_ref->{$trx_id}};
        # Return an error code if transcript has only one exon and --include-monoexonic is not set
        return 2 if scalar @trx_entry == 3 && ! $monoexonic;
		
		#---> Find starting point of alignment and first fragment overlap <---#
		for my $idx ( 2 .. (scalar @trx_entry - 1) ) {

			# For the first fragment, the start coordinate of the alignment should be smaller than the cumulative exon length of the current exon			
			if ( $start_trx <= $trx_entry[$idx]->[3]) {
				# Check if the end coordinate of the alignment is also smaller than the cumulative length of the current exon
				if ( $end_trx <= $trx_entry[$idx]->[3] ) {
					# Return an error code if whole alignment is contained in one exon and --inlcude-monoexonic is not set
					return 2 unless $monoexonic;
					# Calculate genomic start position for alignments to the Watson ("+") and Crick ("-") strands
					$start_gen = $trx_entry[1] eq "+" ? $trx_entry[$idx]->[1] - ($trx_entry[$idx]->[3] - $start_trx) : $trx_entry[$idx]->[1] + ($trx_entry[$idx]->[3] - $end_trx);
				}
				else {
					# Calculate overlap
					my $start_overlap = $trx_entry[$idx]->[3] - $start_trx + 1;
					# Return an error code if overlap is smaller than allowed
					return 3 if $start_overlap < $min_overlap;
					# Push overlap to fragments array
					push @frags, $start_overlap;		
					# Calculate genomic start position for alignments to the Watson ("+") strand
					$start_gen = $trx_entry[$idx]->[1] - $start_overlap + 1 if $trx_entry[1] eq "+";
                                        # Set skip switch
                                        $skip = 0;
				}				
				# Set previous cumulative length for next iteration
				$prev_cum_length = $trx_entry[$idx]->[3];
				# Set flag to indicate that the start coordinate was already found
				$next_index = $idx + 1;
				# Break out of loop
				last;
			}
		}

		#---> Find intermediate and last fragments <---#
		unless ($skip) {
			for my $idx ( $next_index .. (scalar @trx_entry - 1) ) {
	
				# For the last fragment, the end coordinate of the alignment should be smaller than the cumulative exon length of the current exon
				if ( $end_trx <= $trx_entry[$idx]->[3] ) {
					# Calculate overlap
					my $end_overlap = $end_trx - $prev_cum_length;
					# Return an error code if overlap is smaller than allowed
					return 3 if $end_overlap < $min_overlap;
					# Push exon overlap and spanned intron length to fragments array
					push @frags, $trx_entry[$idx]->[4], $end_overlap;
					# Calculate genomic start position for alignments to the Crick ("-") strand
					$start_gen = $trx_entry[$idx]->[0] - $end_overlap + 1 if $trx_entry[1] eq "-";
					# Break out of loop
					last;
				}	
				else {
					# Save length of spanned intron and exon overlap in fragments array
					push @frags, $trx_entry[$idx]->[4], $trx_entry[$idx]->[2];				
					# Set previous cumulative length for next iteration
					$prev_cum_length = $trx_entry[$idx]->[3];				
				}
			}
		}

		#---> Throw warning and return an error code if coordinates are out of bounds (incomplete/incompatible lookup table) <---#
		warn "[WARNING] Alignment (start: $start_trx, end: $end_trx) out of bounds in transcript '$trx_id'! Check integrity and compatibility of exon/transcript annotations used here and during mapping. Entry skipped.\n" and return 1 unless $start_gen;
	
		#---> Reverse fragment order for features on the Crick ("-") strand <---#
		@frags = reverse @frags if $trx_entry[1] eq "-";
	
		#---> Build return array <---#
		push @return, @trx_entry[0..1], $start_gen, @frags;
			
	#---> RETURN VALUE <---#
	return \@return;
	
}
#-----------------------#
sub complement {
### Function: Returns the complement of the input sequence
### Accepts: 1. String (all characters but A, C, G, T and their lower case versions are ignored)
### Returns: Complement of input sequence
### Dependencies: n/a
### Type: Generic

	#---> PASS ARGUMENTS ---#	
	my $seq = shift;

	#---> BODY <---#

		# Build complement by transliteration
		$seq =~ tr/aAcCgGtT/tTgGcCaA/;

	#---> RETURN VALUE <---#
	return $seq;

}
#-----------------------#
sub reverse_CIGAR {
### Function: Reverses a CIGAR string
### Accepts: 1. CIGAR string
### Returns: Reversed CIGAR string
### Dependencies: n/a
### Type: Generic

	#---> PASS ARGUMENTS ---#	
	my $CIGAR = shift;

	#---> SUBROUTINE VARIABLES <---#
	my @CIGAR_new;		

	#---> BODY <---#

		#---> Transform CIGAR string to array <---#
		my @CIGAR_old = split /(\D)/, $CIGAR;

		#---> Pairwise reversal <---#
		push @CIGAR_new, (splice @CIGAR_old, -2) while @CIGAR_old;
		
		#---> Assemble updated CIGAR string <---#
		$CIGAR = join "", @CIGAR_new;		

	#---> RETURN VALUE <---#
	return $CIGAR;

}
#-----------------------#
sub reverse_complement_MD {
### Function: Reverses MD tag
### Accepts: 1. MD tag
### Returns: Reversed MD tag
### Dependencies: n/a
### Type: Generic

	#---> PASS ARGUMENTS ---#	
	my $MD = shift;

	#---> BODY <---#

		#---> Remove tag name and type <---#
		$MD =~ s/MD:Z://;

		# Build complement by transliteration
		$MD =~ tr/aAcCgGtT/tTgGcCaA/;

		#---> Transform MD string to array <---#
		my @MD = split /(\D+)/, $MD;

		#---> Reverse string of deleted letters after caret (MD specifications) <---#
		for my $group (@MD) {
			$group = "^" . reverse $group if $group =~ s/\A\^//;
		}
			
		#---> Reversal <---#
		@MD = reverse @MD;

		#---> Assemble updated MD string <---#
		$MD = join "", "MD:Z:", @MD;

	#---> RETURN VALUE <---#
	return $MD;

}
#-----------------------#
sub add_introns_to_CIGAR {
### Function: For split/spliced alignments, includes introns (Ns) in a CIGAR string 
### Accepts: 1. CIGAR string; 2. array of integers, containing the lengths of exon overlaps and, interspersed, the length(s) of the spanned intron(s), e.g.: 25 (length overlap exon 3), 10000 (length intron between exons 3 and 4), 50 (overlap exon 4), 5000 (length intron between exons 4 and 5), 25 (overlap exon 5)
### Returns: Updated CIGAR string
### Dependencies: n/a
### Type: Generic

	#---> PASS ARGUMENTS ---#	
	my $CIGAR = shift;
	my @frags = @_;

	#---> SUBROUTINE VARIABLES <---#
	my $intron_position = 0;

	#---> BODY <---#

		#---> Transform CIGAR string to array <---#
		my @CIGAR = split /(\D)/, $CIGAR;

		#---> Insert corresponding insert for each pair of fragments <---#		
		for (my $idx_frags = 0; $idx_frags < scalar (@frags - 1); $idx_frags += 2) {
			$intron_position += $frags[$idx_frags];

			#---> Generate array of cumulative lengths that corresponds to CIGAR array; only (mis)matches and deletions count <---#
			my @cum_len = ();	
			my $curr_cum_len = 0;
			for (my $idx = 1; $idx < scalar @CIGAR; $idx += 2) {
				if ($CIGAR[$idx] =~ /^[MD]$/) {	
					$curr_cum_len += $CIGAR[$idx - 1];
				}
				push @cum_len, $curr_cum_len;
			}

			#---> Find insertion position, split previous entry and insert N's <---#
			for (my $idx_cum_len = 0; $idx_cum_len < scalar @cum_len; $idx_cum_len++) {
				if ($intron_position <= $cum_len[$idx_cum_len]) {
					my $part_2 = $cum_len[$idx_cum_len] - $intron_position;
					my $part_1 = $CIGAR[$idx_cum_len * 2] - $part_2;
					if ($part_2) {
						splice @CIGAR, ($idx_cum_len * 2), 2, $part_1, $CIGAR[$idx_cum_len * 2 + 1], $frags[$idx_frags + 1], "N", $part_2, $CIGAR[$idx_cum_len * 2 + 1];
					}
					else {
						splice @CIGAR, ($idx_cum_len * 2), 2, $part_1, $CIGAR[$idx_cum_len * 2 + 1], $frags[$idx_frags + 1], "N";
					}
					last;		
				} 
			}
		}
		
		#---> Assemble updates CIGAR string <---#
		$CIGAR = join "", @CIGAR;

	#---> RETURN VALUE <---#
	return $CIGAR;

}
#-----------------------#
sub rc_bit {
### Function: Extract the 0x10 (sequence reverse complemented) bit of the SAM FLAG field 
### Accepts: 1. SAM FLAG field value
### Returns: Value of 0x10 bit
### Dependencies: n/a
### Type: Generic

	#---> PASS ARGUMENTS ---#	
	my $number = shift;
	
	#---> SUBROUTINE VARIABLES <---#
	my $bin;
	my $bit;

	#---> BODY <---#

		#---> Convert FLAG from decimal to binary <---#
		$bin = sprintf "%012b", $number;

		#---> Insert corresponding insert for each pair of fragments <---#		
		$bit = substr $bin, -5, 1;

	#---> RETURN VALUE <---#
	return $bit;

}

#===================================#
#    DEBUGGING & UNUSED ROUTINES    #
#===================================#
sub print_trx {
	
	my $exon_hoaoa_ref = shift;
	my $id = shift;
	
	my @out_arr = @{$$exon_hoaoa_ref{$id}};		
	print "\nid: " . $id . "\n";
	print "chromosome: " . $out_arr[0] . "\n";
	print "strand: " . $out_arr[1] . "\n";
	for my $idx (2 .. scalar @out_arr - 1) {
		my @in_arr = @{$out_arr[$idx]};
		print "---\nexon: " . ($idx - 2) . "\n";
		print "genome start position: " . $in_arr[0] . "\n";
		print "genome end position: " . $in_arr[1] . "\n";
		print "exon length: " . $in_arr[2] . "\n";
		print "cumulative exon length: " . $in_arr[3] . "\n";
		print "previous intron length: " . $in_arr[4] . "\n";
	}
	print "\n";

}
#-----------------------#
sub CIGAR_aligned_unaligned_stretches {

	my $CIGAR = shift;

	my @CIGAR = split /(\D)/, $CIGAR;

	my @results;
	
	my $last_index = 0;
	
	my( @indices )= grep { $CIGAR[$_] eq "N" } 0..$#CIGAR;

	foreach my $index (@indices) {
		my @temp_array = ();
		@temp_array = @CIGAR[$last_index .. ($index - 2)];
		$last_index = $index + 1;
		my $cum_len = 0;
		for (my $idx = 1; $idx < scalar @temp_array; $idx += 2) {
			if ($temp_array[$idx] =~ /^[MD]$/) {	
				$cum_len += $temp_array[$idx - 1];
			}
		}
		push @results, $cum_len, $CIGAR[$index - 1];
	}

	my @temp_array = ();
	@temp_array = @CIGAR[$last_index .. $#CIGAR];
	my $cum_len = 0;
	for (my $idx = 1; $idx < scalar @temp_array; $idx += 2) {
		if ($temp_array[$idx] =~ /^[MD]$/) {	
			$cum_len += $temp_array[$idx - 1];
		}
	}
	push @results, $cum_len;

	return \@results;

}
#-----------------------#
sub exons_hoaoa_remove_single_exons {
### Function: Removes transcripts with only a single exon entry from exons hoaoa generated by sub "exon_bed_to_hoaoa"
### Accepts: Reference to hash of arrays of arrays
### Returns: Reference to hash of arrays of arrays
### Dependencies: Subroutine "exon_bed_to_hoaoa" (Alexander Kanitz)
### Type: Specialized
	#---> PASS ARGUMENTS ---#
	my $hoaoa_ref = shift;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Removing transcripts with only one exon..." . "\n" unless $quiet;

	#---> BODY <---#

		#---> Remove transcripts with only one exon <---#
		## Traverse through each key $outer_array of hash %$hoaoa_ref
		foreach my $transcript ( keys %$hoaoa_ref ) {
			## Delete hash entry if outer array has only three elements (i.e. one inner array ref / one exon)
			unless (scalar @{$hoaoa_ref->{$transcript}} > 3) {
				delete $hoaoa_ref->{$transcript};
			}
		}
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Single exon transcripts removed." . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return $hoaoa_ref;
	
}
#=======================#
#    SUBROUTINES END    #
#=======================#
