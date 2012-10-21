#!/usr/bin/perl

my @taxa = qw(100 1000 10000 50000 100000);
my @assem = qw(1 10 100 300 500 1000 10000);
my @reads = qw(10 100 1000 10000 100000 200000);

my $j = 5;
my $total = @taxa * @assem * @reads * $j;
my $count = 0;

die unless open(OUT, ">../autovalidate_output/raw.csv");
print OUT "number of taxa\tsize of assemblage\tnumber of reads\tfalse positive\tfalse negative\tminspec false negative\tfalse taxa removed";
foreach my $taxnum (@taxa) {
  foreach my $assemnum (@assem) {
    foreach my $readnum (@reads) {
      my $i;
      for ($i = 1; $i <= $j; ++$i) {
        ++$count;
        print STDERR "\n[$count of $total]\t$taxnum taxon, $assemnum assemblage, $readnum reads, repeat #$i";
        if ($assemnum > $taxnum) {
          print STDERR " - SKIPPING for assem < taxa";
          next;
        }
        my $return = `perl validate_minspec.pl -t $taxnum -a $assemnum -r $readnum 2>&1`;
        my @return = &parse_output($return);
        my $return = join("\t", @return);
        print OUT "\n" . $taxnum . "\t" . $assemnum . "\t" . $readnum . "\t" . $return;
      }
    }
  }
}
close OUT;

sub parse_output {
  my $output = $_[0];
  die $output unless $output =~ /False positive rate: ([^\n]+)%/;
  my $fp = $1;
  die unless $output =~ /False negative rate \(including false negatives due to under-sampling of rare taxa\): ([^\n]+)%/;
  my $fn = $1;
  die unless $output =~ /False negative rate \(only including false negatives generated by minspec\): ([^\n]+)%/;
  my $mfn = $1;
  die unless $output =~ /Correctly removed (\d+) of (\d+)/;
  my $ftr = $1 / $2;
  return ($fp, $fn, $mfn, $ftr);
}
