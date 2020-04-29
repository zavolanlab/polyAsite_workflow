use strict;
use warnings;
use Getopt::Long;

my $help;
my $debug;
my $strict;
my $exclude;
my $min_align  = 4;
my $correction = 0;
GetOptions(
  "strict"       => \$strict,
  "exclude=s"    => \$exclude,
  "correction=i" => \$correction,
  "debug"        => \$debug,
  "min_align=i"  => \$min_align,
  "help"         => \$help,
  "h"            => \$help
);

my $showhelp = 0;
$showhelp = 1 if ( not defined $ARGV[0] );

$showhelp = 1 if ( defined $help );

if ($showhelp) {
  print "[ERROR] Usage: perl $0 --exclude=M:Y input.bed.gz\n\n";
  exit(2);
}

if ( $ARGV[0] !~ /gz$/ ) {
  print STDERR "[ERROR] The input $ARGV[0] is required to be gzipped\n\n";
  exit(2);
}

my %exclude_chr = ();
if ( defined $exclude) {
  my @ex_ar = split( ":", $exclude );
  foreach (@ex_ar) {
    $exclude_chr{$_} = 1;
  }
}

# %ordered stores an ordered set of chromosome names
# starting with chr1 ... chrM
my %ordered     = ();
my $idx_ordered = 0;
foreach my $strand ( '+', '-' ) {
  foreach my $i ( 1 .. 100, 'X', 'Y' ) {
    $idx_ordered++;
    $ordered{"$i:$strand"} = $idx_ordered;
    $idx_ordered++;
  }
}

my %hash      = ();
my %hashdebug = ();
open( F, "gzip -dc $ARGV[0] |" ) or die "Can't open pipe to $ARGV[0]\n";
while (<F>) {
  chomp;
  my @F         = split(/\t/);
  my $readname  = $F[3]; # this is always 1
  my $readcount = 1; # each line corresponds to a distinct read, readcounts not specified expicitly
  next if ( defined $exclude_chr{ $F[0] } );
  if ($strict) {
    my $alignflag = 0;
    if ( $F[5] eq '+' ) {
      if ( $F[6] =~ m/^\d+$/ ) {
        # perfect match: okay
        $alignflag = 1;
      } elsif ( $F[6] !~ m/\d+$/ ) {
        # end badly aligned: skip
      } else {
        if ( $F[6] =~ m/(.*\D)(\d+)$/ ) {
          if ( $2 >= $min_align ) {
            $alignflag = 1;
          }
        }
      }
    } else {
      if ( $F[6] =~ m/^\d+$/ ) {
        # perfect match: okay
        $alignflag = 1;
      } elsif ( $F[6] !~ m/^\d+/ ) {
        # end badly aligned: skip
      } else {
        if ( $F[6] =~ m/^(\d+)(\D.*)/ ) {
          if ( $1 >= $min_align ) {
            $alignflag = 1;
          }
        }
      }
    }
    next if ( $alignflag == 0 );
  }

  my $key = "$F[0]:$F[5]";

  # necessary for the minus strand
  $F[1]++;

  # add the chromosome/strand key if not
  # already in %ordered
  if ( not defined $ordered{$key} ) {
    $idx_ordered++;
    $ordered{$key} = $idx_ordered;
  }

  $hash{$key} = {} if ( not defined $hash{$key} );
  if ($debug) {
    $hashdebug{$key} = {} if ( not defined $hashdebug{$key} );
  }
  if ( $F[5] eq '+' ) {
    $hash{$key}->{ $F[2] } += $readcount;
    if ($debug) {
      $hashdebug{$key}->{ $F[2] } = [] if ( not defined $hashdebug{$key}->{ $F[2] } );
      push @{ $hashdebug{$key}->{ $F[2] } }, $F[3];
    }
  } else {
    $hash{$key}->{ $F[1] } += $readcount;
    if ($debug) {
      $hashdebug{$key}->{ $F[1] } = [] if ( not defined $hashdebug{$key}->{ $F[1] } );
      push @{ $hashdebug{$key}->{ $F[1] } }, $F[3];
    }
  }
}
close(F);

my $idx = 0;
foreach my $key ( sort { $ordered{$a} <=> $ordered{$b} } keys %ordered ) {
  next if ( not defined $hash{$key} );
  my ( $chr, $strand ) = split( /:/, $key );
  foreach my $site ( sort { $a <=> $b } keys %{ $hash{$key} } ) {
    my $start = $site + $correction - 1;
    my $end   = $site + $correction;
    if ( $strand eq '-' ) {
      $start = $site - $correction - 1;
      $end   = $site - $correction;
    }
    my $score = $hash{$key}->{$site};
    $idx++;
    my $additional = '';
    if ($debug) {
      $additional = "\t" . join( ",", @{ $hashdebug{$key}->{$site} } );
    }
    print "$chr\t$start\t$end\t$idx\t$score\t$strand$additional\n";
  }
}

