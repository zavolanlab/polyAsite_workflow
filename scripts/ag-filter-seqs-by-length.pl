use strict;
use warnings;
use Getopt::Long;

my $help;
my $max;
my $min = 1;

GetOptions(
  "max=i" => \$max,
  "min=i" => \$min,
  "help"  => \$help,
  "h"     => \$help
);

my $showhelp = 0;
$showhelp = 1 if ( not defined $max );

if ($showhelp) {
  print STDERR "Usage $0 --max 46 --min 1\n\n";
  exit;
}

my $header = '';
while (<STDIN>) {
  chomp;
  if ( $_ =~ m/^>/ ) {
    $header = $_;
  } else {
    if ( length($_) <= $max ) {
      if ( $header eq '' ) {
        print STDERR "[ERROR] no header found.\n";
        exit;
      }
      print "$header\n$_\n";
      $header = '';
    }
  }
}
