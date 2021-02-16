# mirnatarget

This repository contains scripts that parse the output of miRNA sequences aligned with ssearch against reference bacterial genomes. ssearch can be found in the fasta package (https://github.com/wrpearson/fasta36). The ssearch command used to align miRNA sequences is the following:

For forward strand alignment:
```
ssearch36 -f -16 -g -10.0 -E 3 -T 8 -b=200 -d=200 -r +15/-10 -n -U -W 10 <rev_comped_mirna_sequences> <reference_db> > <ssearch_output>
```

For reverse strand alignment:
```
ssearch36 -f -16 -g -10.0 -E 3 -T 8 -b=200 -d=200 -r +15/-10 -n -U -W 10 <mirna_sequences> <reference_db> > <ssearch_output>
```

Then, parse the alignments using the procedure implemented in psRNATarget, itself inspired from Fahlgren & Carrington (2000).
```
./parseMiRNATargets.pl --help
./parseMiRNATargets.pl --infile ssearch_output.txt > ssearch_parsed.tsv
```

