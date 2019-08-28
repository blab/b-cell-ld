# Analysis of linkage disequilibrium for Dale et al 2019

These scripts and data produce Figure 2 K-P.

Citation: Dale GA, Wilkins DJ, Bohannon CD, Dilernia D, Hunter E, Bedford T, Antia R, Sanz I, Jacob
J. 2019. Clustered mutations at the murine and human IgH locus exhibit significant linkage
consistent with templated mutagenesis. J Immunol: ji1801615.

# Generate LD statistics

Generate R^2 LD statistics for each dataset by running the `ld-analysis.py` script as:

```
python2 ld-analysis.py --dataset rPA_day_8
python2 ld-analysis.py --dataset rPA_day_16
python2 ld-analysis.py --dataset rHA_day_16
python2 ld-analysis.py --dataset tas2016
python2 ld-analysis.py --dataset rabbit
python2 ld-analysis.py --dataset chicken
```

This takes FASTA alignments from the `fastas/` directory and populates TSV files into the
`rsq_tables/` directory. These have the headers `isotype`, `gene`, `distance`, `rsq`, where each row
represents a specific allelic comparison. I've versioned `rsq_tables/rPA_day_8.tsv`, etc... for
convenience.

This script is currently Python 2 compatible only.

# Plot scatterplots

These TSV files are plotted using the supplied `ld-plotting.nb` Mathematica notebook resulting in
`figures/ld_plots.png`.
