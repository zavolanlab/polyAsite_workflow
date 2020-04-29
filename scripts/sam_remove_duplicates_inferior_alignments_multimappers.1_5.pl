#!/usr/bin/env perl
use lib "/scicore/home/zavolan/clipz/newClipz7/lib/perl";
#==================#
#   HEADER START   #
#==================#
### Name: sam_remove_duplicates_inferior_alignments_multimappers.pl
### Created: Aug 29, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Requirements: GetOpt::Long
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
my $head = '';
my $in = '';
my $out = '';
my $new_header = '';
my $multi = -1;
my $mm = '';
my $hm = '';
my $options_result = GetOptions (
	'usage|help' => \$usage,
	'quiet' => \$quiet,
	'print-header' => \$head,
	'new-header=s' => \$new_header,
	'keep-mm:i' => \$multi,
	'mm=s' => \$mm,
	'heavy-mm=s' => \$hm,
	#-----------------------#
	'in=s' => \$in,
	'out=s' => \$out
);

## Die if command line parsing was not successful or required arguments are missing
die $usage_info if $usage || !$options_result;
die $usage_info if !$in || !$out; 

## Die if indicated files do not exist
die "[ERROR] File '$in' not found.\n$usage_info" unless -e $in;
die "[ERROR] File '$new_header' not found.\n$usage_info" unless $new_header eq "" || -e $new_header;

# Unset $head switch if $new_header is set
$head = 0 if $new_header;

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print "Starting '$0'...\n" unless $quiet;

#---> BODY <---#
&filter_sam($in, $out, $multi, $mm, $head, $new_header);

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
'Usage: perl ./sam_remove_duplicates_inferior_alignments_multimappers.pl [OPTIONS] --in [SAM] --out [SAM]

Description: From a sorted SAM file, first removes duplicate records (defined by identical entries for the fields QNAME, FLAG, RNAME, POS & CIGAR), then all QNAME duplicates except for the one(s) with the shortest edit distance. Finally, unless the --keep-mm is set, all alignments of queries with the same edit distance, but different coordinates ("multimappers") are discarded.

Arguments:
--in [SAM]	Input SAM file sorted by QNAME [Required]
--out [SAM]	Output SAM file [Required]
--print-header	Print header (keep input file header if --new-header not specified)
--new-header [FILE]	Uses file indicated in argument as header
--keep-mm [INT]	Keep queries with up to INT different alignments. Set INT to "0" to keep all alignments for each query. By default, all alignments of "multimappers" are removed.
--mm [TAB]	Print the QNAMEs and mapping counts of all "multimappers" to TAB (format: QNAME /TAB/ number of mappings; one entry per line).
--heavy-mm [TAB]        Like --mm, but only prints alignments for reads that map more than --keep-mm times.
--usage|help	Show this information and die
--quiet	Shut up!

Notes:
Script requires NM tags (i.e. edit distances) to be present in all records of the input SAM file.
CAUTION: Only marginal validation of the input file type/format performed!

Version 1.5 (2019-03-27)
Written by Alexander Kanitz on 2013-08-29
';
}
#-----------------------#
sub filter_sam {
### Function: From a sorted SAM file, first removes duplicate records (i.e. same entry name and coordinates, specifically the fields: QNAME, FLAG, RNAME, POS & CIGAR), then all QNAME duplicates except for the one(s) with the shortest edit distance, then (optionally) all alignments of "multimappers" (same QNAME, same edit distance, but different coordinates).
### Accepts: 1. Input file [FILE|SAM]; 2. Output file [FILE|SAM]; 3. Multimapper switch: 0 = remove multimappers, 1 = keep multimappers; 4. Output file for multimapper IDs/QNAMEs; 5. Header switch: FALSE = do not print header, TRUE = print header; 6. Header file (prepends SAM records in output)
### Returns: n/a
### Dependencies: n/a
### Type: Specialized
	#---> PASS ARGUMENTS ---#
	my ($in, $out, $multi, $mm, $head, $new_header) = @_;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Filtering SAM file '$in'..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my $regex_header = '^\@[A-Za-z][A-Za-z](\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$';
	my $regex_comment = '/^\@CO\t.*/';
	my $last_line;																										# holds the last line
	my $last_id;																										# holds the last distinct QNAME/read ID
	my @AoH;																											# holds references to hashes containing the field values of the last lines that share the same QNAME/read ID
	my @field_keys = qw/QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL/;										# holds the bareword keys for the hashes containing field values
	my $final_record = 0;
	
	#---> BODY <---#
	
		#---> Open input and output filehandles <---#
		open IN,  "<", $in  or die "[ERROR] Could not open file '$in'!\n";
		open OUT, ">", $out or die "[ERROR] Could not open file '$out'!\n";
		open MM,  ">", $mm  or die "[ERROR] Could not open file '$mm'!\n" if $mm;
		open HM,  ">", $hm  or die "[ERROR] Could not open file '$hm'!\n" if $hm;
		
		#---> Process header lines (assumed at the top of the file!) <---# 
		while (<IN>) {
			if ( /$regex_header/ || /$regex_comment/ ) {
				print OUT if $head;
			}
			else {	
				$last_line = $_;	
				last;			
			}
		}
	
		#---> Reset filehandle position  <---#
		{
			use bytes;
			seek IN, -length($last_line), 1;
		}
	
		#---> Print header from separate file (if indicated in command line)  <---#	
		if ($new_header) {
			open HEAD, "<", $new_header or die "[ERROR] Could not open file 'new_header'!\n";
				while (<HEAD>) {
					print OUT;				
				}				
			close HEAD;			
		}
		
		#---> Traverse non-header lines <---#
		while (my $line = <IN>) {

			#---> Field values to hash <---#
			chomp $line;
			my @field_values = split "\t", $line;
			die "[ERROR] Input file does not look like a valid SAM file!" unless scalar @field_values >= 11;	# Assert presence of at least 11 fields
			my %fields;
			@fields{@field_keys} = @field_values[0 .. 10];
			foreach (@field_values[11 .. $#field_values]) {
				my ($tag, $value) = split ":", $_, 2;
				$fields{$tag} = $value;
			}
			die "[ERROR] Edit distance ('NM tag') missing!" unless defined $fields{"NM"};						# Assert presence of NM tag

			#---> Manage AoH: Grow if QNAMEs identical, else compare AoH entries, print record(s) and reset AoH <---#
			if ( defined $last_id && $fields{"QNAME"} ne $last_id ) {
				die "[ERROR] SAM file appears to be corrupt!" unless scalar @AoH > 0;							# Assert integrity of SAM records
				@AoH = @{&sam_AoH_filter_records_w_ident_QNAME(\@AoH)} if scalar @AoH > 1;
				my @out_lines = @{&sam_AoH_join_records(\@AoH)};
				if ( $multi < 0 ) {
					print OUT $out_lines[0] if scalar @out_lines == 1;											# Print unique mapper entries
				}
				elsif ( $multi == 0 || $multi >= scalar @out_lines ) {															
					print OUT foreach @out_lines;																# Print all or $multi multimappers if requested
				}
				else
				{
				    print HM $last_id . "\t" . scalar @out_lines . "\n";
				}
				print MM $last_id . "\t" . scalar @out_lines . "\n" if $mm && scalar @out_lines > 1;			# Print multimapper QNAMEs if requested
				@AoH = ();
			}
			push @AoH, \%fields;
			$last_id = $fields{"QNAME"};
			$last_line = $line;
			
		}
		
		#---> Account for final record(s) separately due to EOF <---#
		die "[ERROR] SAM file appears to be corrupt or empty!" unless scalar @AoH > 0;						# Assert integrity of SAM records
		@AoH = @{&sam_AoH_filter_records_w_ident_QNAME(\@AoH)} if scalar @AoH > 1;
		my @out_lines = @{&sam_AoH_join_records(\@AoH)};
		if ( $multi < 0 ) {
			print OUT $out_lines[0] if scalar @out_lines == 1;												# Print unique mapper entries
		}
		elsif ( $multi == 0 || $multi >= scalar @out_lines ) {												
			print OUT foreach @out_lines;																	# Print all or $multi multimappers if requested
		}
		else
        {
        	print HM $last_id . "\t" . scalar @out_lines . "\n";
        }
		print MM $last_id . "\t" . scalar @out_lines . "\n" if $mm && scalar @out_lines > 1;				# Print multimapper QNAMEs if requested
	
		#---> Close input and output filehandles <---#
		close OUT;
		close IN;
		close MM if $mm;
		close HM if $hm;

	#---> END BODY <---#
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Written filtered SAM file to '$out'.\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return 0;
}
#-----------------------#
sub sam_AoH_filter_records_w_ident_QNAME {
### Function: Compares records of a SAM file that share the same QNAME/read ID: True duplicates (i.e. QNAME and coordinates are equal) are discarded first. Then all records but the ones with the lowest edit distances are discarded. The remaining entry or entries are returned in an array of hashes (array length of > 1 if read is a "multimapper").
### Accepts: Array of hashes of SAM records with identical QNAME field (as generated by the subroutine 'filter_sam', written by Alexander Kanitz, 29-AUG-2013)
### Returns: Reference to array of hashes
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $AoH_ref = shift;	
	
	#---> BODY <---#

		#---> Remove "true duplicates" (same QNAME, FLAG, RNAME, POS & CIGAR) / KEEP ONE!!! <---#
		my %true_dup;
		for my $hash_ref (@$AoH_ref) {
			my $id_coord = join "", $hash_ref->{"QNAME"}, $hash_ref->{"FLAG"}, $hash_ref->{"RNAME"}, $hash_ref->{"POS"}, $hash_ref->{"CIGAR"};	# make string to unambiguously define mapping position
			$true_dup{$id_coord} = $hash_ref;																									# Add to hash
		}
		$AoH_ref = [values %true_dup];	

		#---> Keep only ones with shortest edit distance <---#	
		my $min_edit_distance;
		my @edit_distances;
		foreach my $hash_ref (@$AoH_ref) {
			my $edit_distance = ($hash_ref->{"NM"} =~ /^.*:(\d+)/)[0];
			$min_edit_distance = $edit_distance if ! defined $min_edit_distance || $edit_distance < $min_edit_distance;
			push @edit_distances, $edit_distance;
		}	
		for ( my $index = $#edit_distances; $index >= 0; $index-- ) {
			splice @$AoH_ref, $index, 1 if $edit_distances[$index] != $min_edit_distance; 
		}		
	
	#---> RETURN VALUE <---#
	return $AoH_ref;
}
#-----------------------#
sub sam_AoH_join_records {
### Function: Joins the fields of SAM records stored in an array of hashes in the proper order (additional tags in alphanumerical order); re-computes or adds the NG tag field
### Accepts: Array of hashes of SAM records with identical QNAME field (as generated by the subroutine 'filter_sam', written by Alexander Kanitz, 29-AUG-2013)
### Returns: Array of strings, one element for each record, appended with a newline character for easy printing
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $AoH_ref = shift;
	
	#---> SUBROUTINE VARIABLES <---#
	my @field_keys = qw/QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL/;																# holds the bareword keys for the hashes containing field values in the right order
	my @records;	
	
	#---> BODY <---#

		foreach my $record (@$AoH_ref) {
                        my $nh = scalar @$AoH_ref;
                        $record->{NH} = "i:${nh}";
			my @fields_ordered;
			foreach my $key (@field_keys) {
				push @fields_ordered, $record->{$key};
				delete $record->{$key};
			}
			foreach my $extra_field (sort keys %$record) {
				push @fields_ordered, ( $extra_field . ":" . $record->{$extra_field} );
			}
			push @records, ( join( "\t", @fields_ordered) . "\n" );
		}

	#---> RETURN VALUE <---#
	return \@records;
}
#=======================#
#    SUBROUTINES END    #
#=======================#

