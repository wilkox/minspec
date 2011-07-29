#About 

`minspec` is a bioinformatic tool for metagenomic studies. It identifies the most parsimonious set of species which explain the output of a blast search.

A common problem in metagenomics is handling reads which have high sequence identity to more than one species. Because a lot of genomic sequence is highly conserved between species, it's often not clear which species a read originated from. The normal solution is to split the "hit" between all species which match above a certain identity, bit score or E-value threshold, but this is not satisfying for a number of reasons, not least that many of the "identified" species are probably not actually there.

`minspec` solves this by computing the *most parsimonious* set of species needed to explain the observed hits. It does this by framing the computation as a linear programming (LP) problem, which is solved the GNU Linear Programming Kit (GLPK).

The output of `minspec` is _not_ a definitive list of species present in a sample, but rather a result to be interpreted. Although validation tests show `minspec` is generally very reliable when working with simulated metagenomes, it is expected that in real experiments it will report both false positive and false negative results.

`minspec` was inspired by and draws heavily from a similar approach to identifying the minimal set of pathways needed to explain an observed proteome:

[Ye, Y. and Doak, T.G.. A parsimony approach to biological pathway reconstruction/inference for genomes and metagenomes. PLoS Computational Biology 2009, vol 5 num 8](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000465)

#Usage
`minspec` is written in `perl` for a linux or other unix-like environment. It depends on the GNU Linear Programming Toolkit (GLPK) `glpsol`, available from http://www.gnu.org/s/glpk/. If you have trouble getting `minspec` to run, contact the author or (better) modify the script to support your platform and submit a patch.

#Validation
Included in this repository is a script written to validate minspec by simulating a metagenomics experiment. Full details are included in the script (validation/scripts/validate\_minspec.pl).

#License and citation
`minspec` is released into the public domain. To the extent possible under law, all copyright and related or neighboring rights are waived and permission is explicitly and irrevocably granted to copy, modify, adapt, sell and distribute it in any way you choose.

Citation is not necessary but would be appreciated. A manuscript is in preparation, so for now please refer to the project repository at https://github.com/wilkox/minspec/. Patches and contributions are encouraged.

#Author
`minspec` was written by David Wilkins: david@wilkox.org, david.wilkins@unsw.edu.au
