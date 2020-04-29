use strict;
use warnings;
use Getopt::Long;

my @files = ();
my $noip;
my $help;
GetOptions(
  "noip"   => \$noip,
  "help"   => \$help,
  "h"      => \$help
);

# process options
my $showhelp = 0;
@files = @ARGV;
$showhelp = 1 if ( $#files == -1 );
$showhelp = 1 if ($help);

if ($showhelp) {
  print STDERR "Usage: perl $0 --noip file1 file2\n";
  exit 1;
}

# Calculate library size for each sample
my @libsize = ();
foreach my $idx ( 0 .. $#files ) {
  my $sum = 0;
  open( F, "gunzip -c $files[$idx] |" );
  while (<F>) {
    chomp;
    my @F = split(/\t/);
    if ($noip) {
      next if ( $F[3] ne 'OK' );  # IP sites are excluded!
    }
    $sum += $F[4];  # Number of reads that support OK sites are summed
  }
  close(F);
  push @libsize, $sum;
}

my $idx = 2;
foreach my $file (@files) {
  $idx++;
  print "#$idx;$file;$libsize[$idx-3]\n";
}

my %sites       = ();
my %chrs        = ();
my $chr_counter = 0;

foreach my $idx ( 0 .. $#files ) {
  open( F, "gunzip -c $files[$idx] |" );
  while (<F>) {
    chomp;
    my @F = split(/\t/);
    if ($noip) {
      next if ( $F[3] ne 'OK' );
    }
    my $key = "$F[0]\t$F[5]";
    if ( not defined $sites{$key} ) {
      $sites{$key} = {};
      $chr_counter++;
      $chrs{$key} = $chr_counter;
    }
    if ( not defined $sites{$key}->{ $F[2] } ) {
      my @tmp = ();
      foreach (@files) {
        push @tmp, 0;
      }
      $sites{$key}->{ $F[2] } = \@tmp;
    }
    $sites{$key}->{ $F[2] }->[$idx] = $F[4];
  }
  close(F);
}

foreach my $key ( sort { $chrs{$a} <=> $chrs{$b} } keys %chrs ) {
  foreach my $entry ( sort { $a <=> $b } keys %{ $sites{$key} } ) {
    print "$key\t$entry\t", join( "\t", @{ $sites{$key}->{$entry} } ), "\n";
  }
}
