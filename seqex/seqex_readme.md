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
--flank    INT     Length of flanking sequence
--revcomp  STR     (activates reverse complement extraction)
```

`--in` specifies the input FASTA file. The specified sequence must be present in this file.

`--out` specifies the output FASTA file that will contains the sequence region of interest.


## References:

This repository.
