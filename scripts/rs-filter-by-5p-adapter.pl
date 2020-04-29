use strict;
use warnings;
use Getopt::Long;

my $help;
my $adapter;
GetOptions(
  "adapter=s" => \$adapter,
  "help"      => \$help,
  "h"         => \$help
);

my $showhelp = 0;
$showhelp = 1 if ( not defined $adapter );
$showhelp = 1 if ( defined $help );

if ($showhelp) {
  print STDERR "Usage: cat file.fa | $0 --adapter=....TTT\n\n";
  exit 1;
}

$adapter = uc($adapter);
print STDERR "[INFO]   Adapter sequence: $adapter\n";
my $regexp = $adapter;
$regexp =~ s/N/\./g;
print STDERR "[INFO] Regular expression: $regexp\n";

my $header     = '';
my $header_cnt = 0;
my $cnt        = 0;
while (<STDIN>) {
  chomp;
  $cnt++;
  if ( $_ =~ m/^>(.*)/ ) {
    $header = $1;
    $header_cnt = $cnt;
    next;
  }
  if ( $header_cnt + 1 != $cnt ) {
    print STDERR "[ERROR] Script only works with FASTA 1-liners.\n";
    exit;
  }
  if ( $_ =~ m/^($regexp)(.*)/ ) {
    my $newseq = $2;
    if ( length($newseq) > 0 ) {
      print ">$header\n$newseq\n";
    }
  }
}
