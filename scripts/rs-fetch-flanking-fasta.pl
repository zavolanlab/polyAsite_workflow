use strict;
use warnings;
use Getopt::Long;

my $help;
my $genome;
my $upstream   = 100;
my $downstream = 100;
GetOptions(
  "help"         => \$help,
  "h"            => \$help,
  "upstream=i"   => \$upstream,
  "downstream=i" => \$downstream,
  "genome=s"     => \$genome
);

my $showhelp = 0;
$showhelp = 1 if ( defined $help );
$showhelp = 1 if ( not defined $genome );
$showhelp = 1 if ( not defined $ARGV[0] );

if ($showhelp) {
  print STDERR "Usage: perl $0 --genome=genome.fa --upstream=100 --downstream=100 file.bed.gz\n";
  print STDERR "--upstream and --downstream define the region the bed-entry is extended before sequence retrieval\n\n";
  exit 1;
}

if ( not -e $ARGV[0] ) {
  print STDERR "[ERROR] File '$ARGV[0]' not found.\n";
  exit;
}

if ( $ARGV[0] !~ m/.gz/ ) {
  print STDERR "[ERROR] File '$ARGV[0]' does not look like a gzipped file.\n";
  exit;
}

(my $tmp_bed = $ARGV[0]) =~ s/\.gz/\.tmp/;
(my $tmp_fa = $ARGV[0]) =~ s/\.gz/\.tmp\.fa/;

open(OUT, "> $tmp_bed") or die "Can't write to temporary bed-file\n";
open( IN,  "gunzip -c $ARGV[0] |" );
while (<IN>) {
  chomp;
  my @F = split(/\t/);
  my $fa_entry;
  my $start, my $end;
  my $key = "$F[0]:$F[1]:$F[2]:$F[5]";
  if ( $F[5] eq "+" ) {

    $start = $F[1] - $upstream;
    $end = $F[2] + $downstream;
  } else {

    $start = $F[1] - $downstream;
    $end = $F[2] + $upstream;
  }
  print OUT "$F[0]\t$start\t$end\t$key\t1\t$F[5]\n";
}
close(IN);
close(OUT);

`bedtools getfasta -fi $genome -s -bed $tmp_bed > $tmp_fa 2> /dev/null`;

my %seqs   = ();
my $header = '';
open( F, "$tmp_fa" );
while (<F>) {
  chomp;
  if ( $_ =~ m/^>(.+)\:(\d+)\-(\d+)\(([\+-])\)$/ ) {
      my $st;
      my $en;
      if( $4 eq "+" ){
	  $st = $2 + $upstream;
	  $en = $3 - $downstream;
      } else {
	  $st = $2 + $downstream;
	  $en = $3 - $upstream;
      }
    $header = "$1:$st:$en:$4";
  } else {
    if ( $header eq '' ) {
      print STDERR "[ERROR] no header found.\n";
      exit;
    }
    $seqs{ $header } = uc($_);
    $header = '';
  }
}
close(F);


open(BED, "gunzip -c $ARGV[0] |" );
while(<BED>) {
    chomp;
    my @F = split("\t");
    my $key = "$F[0]:$F[1]:$F[2]:$F[5]";
    my $s = ( defined $seqs{ $key } ) ? $seqs{ $key } : "NA";
    print "$_\t$s\n";
}
close(BED);

`rm $tmp_fa`;
`rm $tmp_bed`;



