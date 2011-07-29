#!/usr/bin/perl

#a test to validate minspec

#version 0.1
#written by David Wilkins <david@wilkox.org>, <david.wilkins@unsw.edu.au>
#this software is released into the public domain. To the extent possible under law, all copyright and related or neighboring rights are waived and permission is explicitly granted to copy, modify, adapt, sell and distribute this software in any way you choose.

#this script validates minspec using the following process:
# 1 - randomly generate a set of "taxa", some of which have sequence identity to each other
# 2 - select a subset of these "taxa" to be the "assemblage"
# 3 - simulate a blast of the "assemblage" metagenome against the database of "taxa". to simulate sequence identity between genomes, some "reads" will (at random) have identity to both their true "taxon" and one (or more) "taxa" related to their true "taxon". we will assume that appropriately stringent identity and E-value thresholds have been set on the blast. this blast output will thus represent the common metagenomic situation which minspec is intended to help resolve, where there are high-scoring hits to organisms which are not actually present in the sample, but which have identity to organisms actually in the sample
# 4 - process the "blast output" in minspec
# 5 - compare minspec's minimal species set to the assemblage generated in step 2, and calculate false positive and negative rates

####
## 1 - randomly generate a set of "taxa", some of which have sequence identity to each other

print "\nSTEP 1 - randomly generating 1000 taxa";

#read in lists of adjectives, nouns and animals
my @adjectives = split(/\n/, `cat ../ref/adjectives.list`);
my @nouns = split(/\n/, `cat ../ref/nouns.list`);
my @animals = split(/\n/, `cat ../ref/animals.list`);

#generate database of 1000 "taxa"
print "\nGenerating \"taxa\"...";

until (keys(%taxa) == 1000) {

	my $taxon = @adjectives[int(rand(@adjectives))] . "_" . @nouns[int(rand(@nouns))] . "-" . @animals[int(rand(@animals))];
	next if exists($taxa{$taxon});
	$taxa{$taxon} = "";
	print "\nCreated $taxon";

}

#randomly create sequence identity relationships between taxa
print "\nGenerating relationships between taxa...";

my @taxa = keys(%taxa); #for picking random taxa

foreach $taxon (keys(%taxa)) {

	#max 50 relationships/taxon
	until (keys(%{$relativesOf{$taxon}}) >= 50) {

		#30% of the time, randomly stop adding relationships to the current taxon
		last if rand() > 0.7;

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

		print "\n$taxon and $otherTaxon are now related ($taxon has ". keys(%{$relativesOf{$taxon}}) . " relationships)";
	}

}

####
## 2 - select a subset of these "taxa" to be the "assemblage"

print "\nSTEP 2 - select a subset of 300 taxa to be the assemblage";

#select the subset
