use strict;
use warnings;
use Getopt::Long;

my $help;
my @additional = ();
my $upstream_seq;
my $downstream_seq;
my $total_As;
my $consecutive_As;
my @ds_patterns = ();
GetOptions(
  "help"             => \$help,
  "h"                => \$help,
  "upstream_len=i"   => \$upstream_seq,
  "downstream_len=i" => \$downstream_seq,
  "consecutive_As=i" => \$consecutive_As,
  "total_As=i"       => \$total_As,
  "ds_pattern=s"     => \@ds_patterns,
  "addIPseq:s"       => \@additional
);

my $showhelp = 0;
$showhelp = 1 if ( defined $help );
$showhelp = 1 if ( not defined $ARGV[0] );
$showhelp = 1 if ( not defined $upstream_seq );
$showhelp = 1 if ( not defined $downstream_seq );
$showhelp = 1 if ( not defined $total_As or not defined $consecutive_As );

if ($showhelp) {
  print STDERR
    "Usage: $0 --upstream_seq=0 --downstream_seq=10 --total_As=7 --consecutive_As=6 --ds_pattern=AAAA --ds_pattern=AGAA file.bed.seqs.gz\n";
  print STDERR "\n";
  exit;
}

if ( not -e $ARGV[0] ) {
  print STDERR "[ERROR] File '$ARGV[0]' not found.\n";
  exit 1;
}

my $expected_length = $upstream_seq + $downstream_seq + 1;

open( IN, "gunzip -c $ARGV[0] |" );
while (<IN>) {
  chomp;
  my @F = split(/\t/);

  my $type = 'OK';
  if ( $F[6] eq "NA" or length( $F[6] ) != $expected_length) {
    print STDERR
      "[ERROR] The sequence for the following entry was not retrieved entirely. Entry is marked as IP due to uncertainty.\n$_\n\n";
    $type = "IP";
  } else {
    my $upstream = substr( $F[6], 0,             $upstream_seq );
    my $nuc      = substr( $F[6], $upstream_seq, 1 );
    my $downstream = substr( $F[6], ( $upstream_seq + 1 ), $downstream_seq );

    # either $consecutive_As consecutive Adenosines or at least
    # $total_As Adenosines within the downstream sequence
    my $ip   = 0;
    my $As   = 0;
    my $nt_10_ds = substr( $downstream, 0, 10 );
    foreach my $n ( split( //, $nt_10_ds ) ) {
      $As++ if ( $n eq 'A' );
    }

    foreach my $entry (@additional) {
      $ip = 1 if ( $downstream =~ m/^$entry/ );
    }

    foreach my $entry (@ds_patterns) {
      $ip = 1 if ( $downstream =~ m/^$entry/ );
    }

    $ip = 1 if ( $As >= $total_As );
    $ip = 1 if ( $nt_10_ds =~ m/A{$consecutive_As}/ );
    if ( $ip == 1 ) {
      $type = 'IP';
    }
  }
  $F[3] = $type;
  pop @F;
  print join( "\t", @F ), "\n";
}
