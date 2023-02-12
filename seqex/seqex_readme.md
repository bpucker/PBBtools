# SeqEx: Extraction of sequences from FASTA file

## Usage


```
Usage
python3 seqex3.py --in <FILE> --out <FILE> --contig <STR> --start <INT> --end <INT>

Mandatory:
--in       STR     Input file
--out      STR     Output file
--contig   STR     Sequence name
--start    INT     Start position
--end      INT     End position

Optional:
--flank    INT     Length of flanking sequence[500]
--revcomp  STR     (activates reverse complement extraction)
```

`--in` specifies the input FASTA file. The specified sequence must be present in this file.

`--out` specifies the output FASTA file that will contains the sequence region of interest.

`--contig` specifies the sequence of interest in the input FASTA file.

`--start` specifies the start position of the region of interest on the specified sequence of interest.

`--end` specifies the end position of the region of interest on the specified sequence of interest.

`--flank` specifies the length of flanking sequences that will be extracted in addition to the specified region of interest. Default: 500.

`--revcomp` use of this flag activates the extraction of a sequence as reverse complement. Default: off. 


## References:

This repository.
