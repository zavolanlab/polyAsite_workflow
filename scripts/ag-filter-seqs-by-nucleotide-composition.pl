use strict;
use warnings;
use Getopt::Long;

my $help;
my $max;
my $nuc;

GetOptions(
  "max=f" => \$max,
  "nuc=s" => \$nuc,
  "help"  => \$help,
  "h"     => \$help
);

my $showhelp = 0;
$showhelp = 1 if ( not defined $max );
$showhelp = 1 if ( not defined $nuc );
$showhelp = 1 if ( defined $help );

if ($showhelp) {
  print STDERR "Usage: cat file.fa | $0 --max 2 --nuc N\n\n";
  print STDERR "Usage: cat file.fa | $0 --max 0.8 --nuc A\n\n";
  exit;
}

my $header = '';
while (<STDIN>) {
  chomp;
  if ( $_ =~ m/^>/ ) {
    $header = $_;
  } else {
    my $count = 0;
    if ( $nuc eq 'N' ) {
      $count = ( $_ =~ tr/N// );
    } elsif ( $nuc eq 'A' ) {
      $count = ( $_ =~ tr/A// );
    } elsif ( $nuc eq 'C' ) {
      $count = ( $_ =~ tr/C// );
    } elsif ( $nuc eq 'G' ) {
      $count = ( $_ =~ tr/G// );
    } elsif ( $nuc eq 'T' ) {
      $count = ( $_ =~ tr/T// );
    }
    if ( $max < 1 or $max eq '1.0' ) {
      $count = $count / length($_);
    }
    if ( $count <= $max ) {
      if ( $header eq '' ) {
        print STDERR "[ERROR] no header found.\n";
        exit;
      }
      print "$header\n$_\n";
      $header = '';
    }
  }
}
