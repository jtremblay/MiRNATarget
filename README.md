# MiRNATarget (aka psRNATarget_command_line)

This repository contains an implementation of the psRNATarget software (2018) (https://doi.org/10.1093/nar/gky316) for which there is a web service available here (http://plantgrn.noble.org/psRNATarget/). Using this web service is relevant for small datasets, but not for large datasets in the context of repeated usage by a bioinformatic production pipeline that can take full advantage of high performance computing capabilities. Unfortunately, no code is publicly available for that software and as such, no local implementation with custom databases is possible. psRNATarget in its current state (i.e. 2018 release) actually loosely describes an implementation of a 2011 previous psRNATarget release (https://doi.org/10.1093/nar/gkr319) which in turn refers to the Fahlgren & Carrington (2009) book chapter.  

Here, we present code that implement the psRNATarget methodology. It parses the output of miRNA sequences aligned with SSEARCH against reference bacterial genomes. This script can be deployed and used on a HPC system. SEARCH can be found in the fasta package (https://github.com/wrpearson/fasta36). 

The SSEARCH command used to align miRNA sequences is the following command inspired from TargetFinder (https://github.com/carringtonlab/TargetFinder) - look for line 243 of ```targetfinder.pl``` in particular. By trial and error, we've found the following ssearch commands to exactly reproduce the output of psRNATarget.

For forward strand alignment:
```
ssearch36 \
 -f -8 -g -2 -E 10000 -T 8 -b 200 -r +4/-3 -n -U -W 10 -N 20000 \
 <rev_comped_mirna_sequences> \
 <reference_db> \
 > <ssearch_output>
```

For reverse strand alignment:
```
ssearch36 \
-f -8 -g -2 -E 10000 -T 8 -b 200 -r +4/-3 -n -U -W 10 -N 20000 \
 <mirna_sequences> \
 <reference_db> \
 > <ssearch_output>
```
Note that in our tests, the ```-r``` argument had to be absolutely set to ```-r +4/-3```. The ```-f``` and ```-g``` parameters that were found to work are the following:
```-f -8 -g -3```, ```-f -8 -g -2``` or ```-f -9 -g 2```.
Then, simply parse the alignments using the procedure implemented in psRNATarget, itself inspired from Fahlgren & Carrington (2000):
```
./parseMiRNATargets.pl --help
./parseMiRNATargets.pl --infile ssearch_output.txt > ssearch_parsed.tsv
```
We compared the results of this script with the ones given by psRNATarget and it gives identical results. 

If you use MiRNATarget in your work, please cite:

Tremblay, Julien

MiRNATarget 1.0 : miRNA target finder

https://github.com/jtremblay/mirnatarget
