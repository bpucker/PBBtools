# Heatmapper: Visualization of gene expression

## Usage

```
Usage
python3 heatmap_plotter.py --genes <FILE> --exp <FILE> --out <FILE>

Mandatory:
--genes    STR   Candidate gene file.
--exp      STR   Count table.
--out      STR   Output folder

Optional:
--norm (activates normalization per gene)
```

`--genes` specifies a text file containing the genes of interest. Each line lists one gene ID. These IDs need to match the IDs in the first column of the count table. The first column can be followed by additional columns with trivial names in a second column. Columns must be separated by tabs.

`--exp` specifies a text file that contains all the expression data (count table, matrix). Gene IDs are given in the first column and sample IDs are given in the first row.

`--out` specifies the heatmap output file. The file extension can determine the file type. File types require support from Python and the operating system. Frequently used file types PNG, JPG, SVG, and PDF.

`--norm` this flag activates a normalizaiton per gene. Default: off.


## References
This repository.
