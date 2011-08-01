#About 

`minspec` is a simple bioinformatic tool for metagenomic studies. It identifies the most parsimonious set of species which explain the output of a blast search.

A common problem in metagenomics is handling reads which have high sequence identity to more than one species. Because a lot of genomic sequence is conserved between species, it's often unclear which species a read originated from. The normal solution is to split hits between all species which match above a certain identity, bit score or E-value threshold, but this is not satisfying for several reasons, not least that many of the identified species are probably not actually there.

`minspec` helps fix this by computing the *most parsimonious* set of species needed to explain the observed hits. It does this by framing the computation as a linear programming (LP) problem, which is solved the GNU Linear Programming Kit (GLPK).

The output of `minspec` is _not_ a definitive list of species present in a sample, but rather a result to be interpreted. Although tests show `minspec` is generally very reliable when working with simulated metagenomes, it is expected that in real experiments it will report both false positive and false negative results.

`minspec` was inspired by and draws heavily from a similar approach to identifying the minimal set of pathways needed to explain an observed proteome:

[Ye, Y. and Doak, T.G.. A parsimony approach to biological pathway reconstruction/inference for genomes and metagenomes. PLoS Computational Biology 2009, vol 5 num 8](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000465)

#Usage
`minspec` is written in `perl` for a linux or other unix-like environment. It depends on the GNU Linear Programming Toolkit (GLPK) `glpsol`, available from http://www.gnu.org/s/glpk/. If you have trouble getting `minspec` to run, contact the author or (better) modify the script to support your platform and submit a patch.

	USAGE:
	
	minspec -b <blast hittable>

	OPTIONAL:

	-g      <filename> Produce a BLAST hit table containing only hits to species present in the minimal set, suitable for processing with GAAS. Default filename is <blast hittable>.filtered
	-l      <filename> Produce a list of species, indicating whether or not they are present in the minimal set ('1' = present, '0' = not present). Default filename is <blast hittable>.minimal.list
	-max    <maximum read count> Set a maximum number of reads with identity to a species, less than which the species may still be parsimoniously eliminated but equal to or more than which the species will be marked as present regardless of its presence in the minimal set. Set to 50 if -m flag is provided without a value.

#Validation
Included in this repository is a script written to validate minspec by simulating a metagenomics experiment. Full details are included in the script (validation/scripts/validate\_minspec.pl).

#License and citation
`minspec` is released into the public domain. To the extent possible under law, all copyright and related or neighboring rights are waived and permission is explicitly and irrevocably granted to copy, modify, adapt, sell and distribute it in any way you choose.

Citation is not necessary but would be appreciated. A manuscript is in preparation, so for now please refer to the project repository at https://github.com/wilkox/minspec/. Code contributions of any kind are encouraged.

#Author
`minspec` was written by David Wilkins: david@wilkox.org, david.wilkins@unsw.edu.au
