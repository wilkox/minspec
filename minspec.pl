#!/usr/bin/perl

#minspec: determines minimal set of species needed to explain 
# species assignments from a dataset of metagenomic reads, 
# eliminating spurious species assignments.

#based on the approach of, and borrows heavily from, MinPath:
# Ye, Y. and Doak, T.G.. A parsimony approach to biological pathway 
# reconstruction/inference for genomes and metagenomes. PLoS
# Computational Biology 2009, vol 5 num 8
#http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000465

#written by David Wilkins <david@wilkox.org>
#minspec lives at https://github.com/wilkox/minspec

#this software is released into the public domain. To the extent 
# possible under law, all copyright and related or neighboring
# rights are waived and permission is explicitly and irrevocably
# granted to copy, modify, adapt, sell and distribute this software
# in any way you choose.

$USAGE = q/USAGE:

minspec -b <blast hittable>

OPTIONAL:

-g	<filename> Produce a BLAST hit table containing only hits to species present in the minimal set, suitable for processing with GAAS. Default filename is <blast hittable>.filtered
-l	<filename> Produce a list of species, indicating whether or not they are present in the minimal set ('1' = present, '0' = not present). Default filename is <blast hittable>.minimal.list
-max	<maximum read count> Set a maximum number of reads with identity to a species, less than which the species may still be parsimoniously eliminated but equal to or more than which the species will be marked as present regardless of its presence in the minimal set. Set to 50 if -m flag is provided without a value.
/;

use Getopt::Long;

GetOptions (
'b=s' => \$blastoutputfile,
'g:s' => \$makegaas,
'l:s' => \$makelist,
'max:s' => \$maxthreshold,
) or die("$USAGE");

#check for required arguments and set defaults

die ("$USAGE\n") if !defined $blastoutputfile;
die ("ERROR - you did not specify any output format, so there's no use in running the script! Try using -g or -l\n") unless (defined $makegaas || defined $makelist);
$maxthreshold = 50  if (!$maxthreshold && defined $maxthreshold);
$makegaas = "$blastoutputfile.filtered" if (!$makegaas && defined $makegaas);
$makelist = "$blastoutputfile.minimal.list" if (!$makelist && defined $makelist);

##BODY
&doLP;
&makegaas if defined $makegaas;
&makelist if defined $makelist;
exit;
##END BODY

#SUBS

sub doLP { #set up and run the lp

	#read in blast output
	die ("ERROR - could not open blast output file $blastoutputfile\n") unless open(BLAST, "<$blastoutputfile");

	#set unique ids
  #the MPS file format has some very strict requirements -
  # replacing read and OTU names with these IDs ensures these
  # are never inadvertantly broken
	$readuid = "AAAAAAA";
	$specuid = "LLLLLL";

  #get read-species mappings from the BLAST output
	while ($line = <BLAST>) {
		chomp $line;
		unless ($line =~ /^(\S+)\s+(\S+)/) {
			print ("ERROR - malformed line in blast output file $blastoutputfile, line $.\n");
			next;
		}
		$read = $1;
		$species = $2;

		#assign unique IDs to read and OTU if
    # they don't already have them
		if (!exists ($readuidof{$read})) {
			$readuidtrans{$readuid} = $read;
			$readuidof{$read} = $readuid;
			}
		unless (exists ($specuidof{$species})) {
			$specuidtrans{$specuid} = $species;
			$specuidof{$species} = $specuid;
			}

    #store read-species mappings
		push(@{$allreads{$readuidof{$read}}}, $specuidof{$species});
		push(@{$allspecies{$specuidof{$species}}}, $readuidof{$read});
		++$readuid;
		++$specuid;

		}
	close BLAST;

	#produce MPS file
	die unless open(MPS, ">test.mps");
	print MPS q/NAME          PATH
ROWS
 N  NUM/;

	#first, list all the reads
	foreach $read (keys(%allreads)) {
		print MPS "\n G  $read";
		}

	#next, list all species with their read mappings
	print MPS "\nCOLUMNS";
	foreach $species (keys(%allspecies)) {
		$speciesspacerlength = 10 - length($species);
		print MPS "\n    $species" . " "x$speciesspacerlength . "NUM                1";
		undef %preventduplicatereads;
		foreach $read (@{$allspecies{$species}}) {
			next if exists($preventduplicatereads{$read});
			$readspacerlength = 19 - length($read);
			print MPS "\n    $species" . " "x$speciesspacerlength . "$read" . " "x$readspacerlength . "1";
			$preventduplicatereads{$read} = "";
			}
		}

	#list all reads (constraint function)
	print MPS "\nRHS";
	foreach $read (keys(%allreads)) {
		$readspacerlength = 17 - length($read);
		print MPS "\n    RHS1      $read" . " "x$readspacerlength . "1.0";
		}

	#list all species (not sure why this is needed by GLPSOL complains
  # without it)
	print MPS "\nBOUNDS";
	foreach $species (keys(%allspecies)) {
		print MPS "\n BV BND1      $species    ";
		}

	#and done
	print MPS "\nENDATA\n";
	close MPS;

	#run the LP
	die ("ERROR - LP execution with glpsol failed (command: glpsol test.mps -o test.mps.LPout)\n") unless system("glpsol test.mps -o test.mps.LPout") == 0;

	#parse the LP output
	die ("ERROR - could not open the LP output test.mps.LPout\n") unless open(LP, "<test.mps.LPout");
	while ($line = <LP>) {
		$reading = 1 if $line =~ /Column name/;
		next unless $reading == 1;
		next if $line =~ /Column name/;
		chomp $line;
		next if $line =~ /^-/;
		last if $line eq "";
		@line = split(/\s+/, $line);
		$species = @line[2];
		$speciesname = $specuidtrans{$species};
		$presence{$speciesname} = @line[4]; #set to 1 if the species is present in the minimal set, 0 if not

		if (defined $maxthreshold && $presence{$2} == 0 && @{$allspecies{$specuidof{$speciesname}}} >= $maxthreshold) { #if a max threshold is set AND the species presence is set to zero AND the count of reads hitting to that species is greater than the max threshold...
			$presence{$speciesname} = 1; #...set species to present
			}

		}

	#clean up
	close LP;
	system("rm test.mps test.mps.LPout");

} #end of doLP sub

sub makegaas { #makes an blast hittable for gaas, identical to the input but with non-parsimonious species removed

	die ("ERROR - could not open blast output file $blastoutputfile\n") unless open(BLAST, "<$blastoutputfile");
	die ("ERROR - could not create filtered blast hittable for gaas at $makegaas\n") unless open(GAAS, ">$makegaas");

	while ($line = <BLAST>) {
		chomp $line;
		die ("ERROR - malformed blast hittable line at line $. of $blastoutputfile\n") unless $line =~ /^(\S+)\s+(\S+)/;
		print GAAS "$line\n" if $presence{$2} == 1; #print the line if this species is present in the minimal set
		}

	close BLAST;
	close GAAS;

} #end of makegaas sub

sub makelist { #makes a simple list of species, indicating whether or not they are present in the minimal set

	die ("ERROR - could not create list of species at $makelist\n") unless open(LIST, ">$makelist");

	foreach $species (sort {$a <=> $b} (keys(%presence))) {
		print LIST "$species\t$presence{$species}\n";
		}

	close LIST;

} #end of makelist sub
