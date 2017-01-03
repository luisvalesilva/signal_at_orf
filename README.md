# signal_at_orf
Calculate signal at ORF genome-wide. This is a python 3 implementation of the namesake function in [hwglabr](https://github.com/hochwagenlab/hwglabr).

Given wiggle data generated using the lab's ChIP-seq analysis pipelines, this function pulls out the ChIP signal over all ORFs in the yeast genome (flanked by half the ORF length on either side) scaled to a constant length.

The [plot_signal_at_orf.py](plot_signal_at_orf.py) script allows a quick visalization of the output. Given one or two TSV files it calculates mean signal grouped by position and plots the result.

#### For command line usage run:
```
python signal_at_orf.py -h
python plot_signal_at_orf.py -h
```

