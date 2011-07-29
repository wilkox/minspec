#!/usr/bin/perl

#a test to validate minspec

#version 0.1
#written by David Wilkins <david@wilkox.org>, <david.wilkins@unsw.edu.au>
#this software is released into the public domain. To the extent possible 
#under law, all copyright and related or neighboring rights are waived and 
#permission is explicitly granted to copy, modify, adapt, sell and 
#distribute this software in any way you choose.

#this script validates minspec using the following process:
# 1 - randomly generate a set of "taxa", some of which have 
#	sequence identity to each other
# 2 - select a subset of these "taxa" to be the "assemblage"
# 3 - simulate a blast of the "assemblage" metagenome against 
#	the database of "taxa". to simulate sequence identity between 
#	genomes, some "reads" will (at random) have identity to both their 
#	true "taxon" and one (or more) "taxa" related to their true "taxon". 
#	we will assume that appropriately stringent identity and E-value 
#	thresholds have been set on the blast. this blast output will thus 
#	represent the common metagenomic situation which minspec is intended 
#	to help resolve, where there are high-scoring hits to organisms 
#	which are not actually present in the sample, but which have 
#	identity to organisms actually in the sample
# 4 - process the "blast output" in minspec
# 5 - compare minspec's minimal species set to the assemblage generated 
#	in step 2, and calculate false positive and negative rates

#verbosity
if (@ARGV[0] eq "--verbose" || @ARGV[0] eq "-v") {
	$verbose = 1;
}
print "\nRun with --verbose or -v for verbose output" unless $verbose == 1;

#make sure working dirs exist
mkdir "../blast_output" unless -d "../blast_output";
mkdir "../minspec_output" unless -d "../minspec_output";

####
## 1 - randomly generate a set of "taxa", some of which have sequence identity to each other

print "\nSTEP 1 - randomly generating 50000 taxa" if $verbose == 1;

#read in lists of adjectives, nouns and animals
my @adjectives = split(/\n/, `cat ../ref/adjectives.list`);
my @nouns = split(/\n/, `cat ../ref/nouns.list`);
my @animals = split(/\n/, `cat ../ref/animals.list`);

#generate database of 50000 "taxa"
print "\nGenerating \"taxa\"..." if $verbose == 1;

until (keys(%taxa) == 50000) {

	my $taxon = @adjectives[int(rand(@adjectives))] . "_" . @nouns[int(rand(@nouns))] . "-" . @animals[int(rand(@animals))];
	next if exists($taxa{$taxon});

	$taxa{$taxon} = "";
	print "\nCreated $taxon" if $verbose == 1;

}

#randomly create sequence identity relationships between taxa
print "\nGenerating relationships between taxa..." if $verbose == 1;

my @taxa = keys(%taxa); #for picking random taxa

foreach $taxon (keys(%taxa)) {

	#max 50 relationships/taxon
	until (keys(%{$relativesOf{$taxon}}) >= 50) {

		#30% of the time, randomly stop adding relationships to the current taxon
		last if rand() > 0.3;

		#pick another taxon at random
		my $otherTaxon = @taxa[int(rand(@taxa))];

		#don't make relationship if other taxon is self, relationship
		#already exists, or other taxon already has 50 relationships
		next if $otherTaxon eq $taxon;
		next if exists $relativesOf{$taxon}{$otherTaxon};
		next if keys(%{$relativesOf{$otherTaxon}}) >= 50;

		#make the relationship (in both directions)
		$relativesOf{$taxon}{$otherTaxon} = "";
		$relativesOf{$otherTaxon}{$taxon} = "";

		print "\n$taxon and $otherTaxon are now related ($taxon has ". keys(%{$relativesOf{$taxon}}) . " relationships)" if $verbose == 1;
	}

}

####
## 2 - select a subset of these "taxa" to be the "assemblage"

print "\nSTEP 2 - select a subset of 300 taxa to be the assemblage" if $verbose == 1;

#select the subset
until (keys(%assemblage) == 300) {
	my $taxon = @taxa[int(rand(@taxa))];
	next if exists $assemblage{$taxon};

	#taxon relative abundance will follow a log curve (roughly)
	my $abundance = 1 / (2.71828183 ** (1 + keys(%assemblage) / 20));
	$assemblage{$taxon} = $abundance;

	print "\n$taxon is part of the assemblage, with a relative abundance of $abundance" if $verbose == 1;
}

####
## 3 - simulate a blast of the "assemblage" metagenome against the database of "taxa". to simulate sequence identity between genomes, some "reads" will (at random) have identity to both their true "taxon" and one (or more) "taxa" related to their true "taxon". we will assume that appropriately stringent identity and E-value thresholds have been set on the blast. this blast output will thus represent the common metagenomic situation which minspec is intended to help resolve, where there are high-scoring hits to organisms which are not actually present in the sample, but which have identity to organisms actually in the sample

print "\nSTEP 3 - simulating blast of a metagenomic sample from the assemblage against the database of taxa\n" if $verbose == 1;

#open the "blast output"
die unless open(OUT, ">../blast_output/validation.blast_output");

#"blast" 100000 "reads"
my @assemblage = keys(%assemblage); #for random picking
my $read = 0; #initial read ID
my $blastline = "100	500	0	0	1	500	1	500	0	500"; #this is the rest of the blast output line - not actually important for minspec but will be added to each line for completeness

until ($read == 100000) {

	print "\rread $read of 100000" if $verbose == 1;

	#randomly select a taxon from the assemblage
	my $taxon = @assemblage[int(rand(@assemblage))];

	#skip taxa at random to simulate abundance curve
	#randomness determined by taxon abundance cutoff
	#set when assemblage generated
	next if rand() > $assemblage{$taxon};

	#with a 30% chance on each iteration, also produce hits to other related taxa.
	#by shifting off one end of the list instead of picking randomly, 
	#we simulate a curve of relationship distance
	my @relatives = keys(%{$relativesOf{$taxon}}); #for random picking
	until (@relatives == 0) {
		last if rand() > 0.3;
		my $otherTaxon = shift(@relatives);
		print OUT "read$read\t$otherTaxon\t$blastline\n";
		$presentInBlast{$otherTaxon};
	}

	#produce a hit to the real taxon - we put
	#this at the end so minspec cannot simply
	#select the best hit (not that it does)
	print OUT "read$read\t$taxon\t$blastline\n";
	$presentInBlast{$taxon};

	++$read;
}

close OUT;

####
## 4 - process the "blast output" in minspec

print "\nSTEP 4 - processing the pseudo-blast output with minspec" if $verbose == 1;

system("perl ../../minspec.pl -b ../blast_output/validation.blast_output -l ../minspec_output/validation_minspec_output.list");

####
## 5 - compare minspec's minimal species set to the assemblage generated in step 2, and calculate false positive and negative rates

print "\nSTEP 6 - comparing minspec-determined minimal set to actual assemblage" if $verbose == 1;

#read in minspec output
my @minspecOutput = split(/\n/, `cat ../minspec_output/validation_minspec_output.list`);

#go though output and compare to assemblage
foreach (@minspecOutput) {
	$_ =~ /^(\S+)\t([1|0])$/;
	my $taxon = $1;
	my $status = $2;

	#for taxa which are part of the minimal set
	if ($status == 1) {
		++$falsePositive unless exists $assemblage{$taxon};
	}

	#for taxa which are not part of the minimal set but
	#were part of the assemblage
	if ($status == 0) {
		++$falseNegativeMinspec if exists $assemblage{$taxon};
	}

	#to calculate the total false negative rate
	$status{$taxon} = $status;
}

#find taxa in assemblage not reported as present by minspec
foreach my $taxon (@assemblage) {
	++$falseNegative unless $status{$taxon} == 1;
}

#report the results!

#false +ve rate is proportion of "false" taxa (present in blast results
#but not in assemblage) reported by minspec as positive
my $falsePositiveRate = abs((100 * $falsePositive) / (keys(%presentInBlast) - @assemblage));

#false -ve rate is proportion of assemblage taxa not
#reported as present by minspec
my $falseNegativeRate = (100 * $falseNegative) / @assemblage;

#the false -ve rate due to minspec is the proportion of taxa
#present in assemblage and blast output but eliminated by minspec
my $falseNegativeMinspecRate = (100 * $falseNegativeMinspec) / @minspecOutput;

print "\nRESULTS:\nFalse positive rate: $falsePositiveRate%\nFalse negative rate (including false negatives due to under-sampling of rare taxa): $falseNegativeRate%\nFalse negative rate (only including false negatives generated by minspec): $falseNegativeMinspecRate%";
