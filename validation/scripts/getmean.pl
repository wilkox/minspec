#!/usr/bin/perl

my %fp;
my %fn;
my %fnm;
my @all;
print "assemblage/taxa\tnumber of reads\tfalse positive\tfalse negative\tminspec false negative";
while (my $line = <STDIN>) {
  next if $. == 1;
  chomp $line;
  die unless $line =~ /^(\d+\t\d+\t\d+)\t([\d|\.]+)\t([\d|\.]+)\t([\d|\.]+)$/;
  push(@all, $1) unless exists $fp{$1};
  $fp{$1} += $2 / 5;
  $fn{$1} += $3 / 5;
  $fnm{$1} += $4 / 5;
}

foreach my $line (@all) {
  die unless $line =~ /^(\d+)\t(\d+)\t(\d+)/;
  my $ratio = $2 / $1;
  print "\n$ratio\t$3\t$fp{$line}\t$fn{$line}\t$fnm{$line}";
}
