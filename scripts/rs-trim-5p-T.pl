use strict;
use warnings;
use Getopt::Long;

my $help;
my $nucs;
my $minLen;

GetOptions(
  "nuc=s"    => \$nucs,
  "help"     => \$help,
  "h"        => \$help,
  "minLen=i" => \$minLen
);

my $nuc = '[' . $nucs . ']';

my $showhelp = 0;
$showhelp = 1 if ( not defined $nuc );
$showhelp = 1 if ( defined $help );
$showhelp = 1 if ( not defined $minLen );

if ($showhelp) {
  print STDERR "Usage: cat sample.fa | perl $0 --nuc=C --minLen=x\n\n";
  exit;
}

my $header = '';
while (<STDIN>) {
  chomp;
  if ( $_ =~ m/^>/ ) {
    $header = $_;
  } else {
    if ( $header eq '' ) {
      print STDERR "[ERROR] no header found.\n";
      exit;
    }
    my $seq = $_;
    $seq =~ s/^$nuc+//g;
    if ( length($seq) >= $minLen ) {
      print "$header\n$seq\n";
    }
    $header = '';
  }
}
