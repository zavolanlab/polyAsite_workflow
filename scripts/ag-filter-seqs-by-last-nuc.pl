use strict;
use warnings;
use Getopt::Long;

my $help;

GetOptions(
  "help" => \$help,
  "h"    => \$help
);

my $showhelp = 0;
$showhelp = 1 if ( defined $help );

if ($showhelp) {
  print STDERR "Usage: zcat file.fasta | $0\n\n";
  exit;
}

my $header = '';
while (<STDIN>) {
  chomp;
  if ( $_ =~ m/^>/ ) {
    $header = $_;
  } else {
    if ( $_ !~ m/A$/ ) {
      if ( $header eq '' ) {
        print STDERR "[ERROR] no header found.\n";
        exit;
      }
      print "$header\n$_\n";
      $header = '';
    }
  }
}

