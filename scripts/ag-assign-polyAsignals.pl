use strict;
use warnings;
use Getopt::Long;

my $help;
my $genome;
my @motif = ();
GetOptions(
  "help"     => \$help,
  "h"        => \$help,
  "genome=s" => \$genome,
  "motif=s"  => \@motif
);

my $showhelp = 0;
$showhelp = 1 if ( defined $help );
$showhelp = 1 if ( not defined $genome );
$showhelp = 1 if ( not defined $ARGV[0] );
$showhelp = 1 if ( $#motif == -1 );

if ($showhelp) {
  print "Usage: perl $0 --genome=genome.fa --motif=AATAAA --motif=ATTAAA file.tsv.gz\n\n";
  exit 1;
}

if ( not -e $ARGV[0] ) {
  print STDERR "[ERROR] File '$ARGV[0]' not found.\n";
  exit;
}

### Added by CJH, 23.10.2018
if ( not -e $genome ) {
  print STDERR "[ERROR] Genome file '$genome' not found.\n";
  exit;
}
###

my $upstream   = 59;
my $downstream = 10;


(my $tmp_bed = $ARGV[0]) =~ s/\.gz/\.tmp/;
(my $tmp_fa = $ARGV[0]) =~ s/\.gz/\.tmp\.fa/;

open(OUT, "> $tmp_bed") or die "Can't write to temporary bed-file\n";
open( IN, "gunzip -c $ARGV[0] |" );
while (<IN>) {
  chomp;
  next if ( $_ =~ m/^#/ );

  my @F = split(/\t/);
  my $fa_entry;
  my $start, my $end;
  my $key = "$F[0]:$F[2]:$F[1]";
  if ( $F[1] eq "+" ) {
    $start = $F[2] - 1 - $upstream;
    $end   = $F[2] + $downstream;
  } else {
    $start = $F[2] -1 - $downstream;
    $end   = $F[2] + $upstream;
  }
  print OUT "$F[0]\t$start\t$end\t$key\t1\t$F[1]\n";
}
close(IN);
close(OUT);

# get valid signals from input
my %motifs = ();
foreach my $sig (@motif) {
  $motifs{$sig} = { 'regexp' => $sig, 'positions' => [], 'length' => 6, 'PAS' => 1 };
}

`bedtools getfasta -fi $genome -s -bed $tmp_bed > $tmp_fa 2> /dev/null`;

my %seqs   = ();
my $header = '';
open( F, "$tmp_fa" );
while (<F>) {
  chomp;
  if ( $_ =~ m/^>(.+)\:(\d+)\-(\d+)\(([\+-])\)$/ ) {
      my $en;
      if( $4 eq "+" ){
	  $en = $3 - $downstream;
      } else {
	  $en = $3 - $upstream;
      }
      $header = "$1:$en:$4";
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

open( IN, "gunzip -c $ARGV[0] |" );
while (<IN>) {
  chomp;
  if ( $_ =~ m/^#/ ) {
    print "$_\n";
    next;
  }
  my @F = split(/\t/);
  my $key = "$F[0]:$F[2]:$F[1]";
  my $sequence = ( defined $seqs{ $key } ) ? $seqs{ $key } : '';

  # reset $motifs{$key}->{positions}
  foreach my $key ( keys %motifs ) {
    $motifs{$key}->{positions} = [];
  }

  # count motif occurences
  for my $key ( keys %motifs ) {
    for my $i ( 0 .. length($sequence) - $motifs{$key}->{length} ) {
      my $submotif = substr $sequence, $i, $motifs{$key}->{length};
      if ( $submotif =~ m/^$motifs{$key}->{regexp}$/i ) {
	  push @{ $motifs{$key}->{positions} }, $i;
      }
    }
  }

  my @string_tmp = ();

  foreach my $key ( keys %motifs ) {
    next if ( $motifs{$key}->{PAS} == 0 );
    next if ( $#{ $motifs{$key}->{positions} } == -1 );

    foreach my $position ( @{ $motifs{$key}->{positions} } ) {

      my $newpos;

      # index balancing necessary:
      # no 0 is possible since first nt upstream of the CS is annotated as -1,
      # first nt downstream of the CS is annotated as +1
      if ( ( $position - 60 ) < 0 ) {
        $newpos = $position - 60;
      } else {
        $newpos = $position - 60 + 1;
      }
      my $realpos = 'NA';
      if ( $F[1] eq '+' ) {
        if ( $newpos < 0 ) {
          $realpos = $F[2] + $newpos + 1;
        } else {
          $realpos = $F[2] + $newpos;
        }
      } else {
        if ( $newpos < 0 ) {
          $realpos = $F[2] - $newpos - 1;
        } else {
          $realpos = $F[2] - $newpos;
        }
      }
      push @string_tmp, "$key\@$newpos@" . $realpos;
    }
  }

  my $PAS = join( ";", sort {get_position($a) <=> get_position($b)} @string_tmp );
  $PAS = 'NA' if ( $PAS eq '' );

  print "$_\t$PAS\t$sequence\n";
}
close(IN);

`rm $tmp_fa`;
`rm $tmp_bed`;

sub get_position {
    my ($v) = @_;
    my $c;
    if($v =~ /\@(\S+)\@/) {
	$c = $1;
    }
    else {
	print STDERR "Misformatted signal $v\n";
	exit(1);
    }
    return $c;
}
