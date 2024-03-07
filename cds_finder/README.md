## Usage ##

This script was developed to discover coding sequences (CDS) within a transcript sequence. For example, this script could be applied to screen a transcriptome assembly. It does not work on a genome sequence, because the CDS must be continuous and introns would prevent the discovery.

```
Usage
python3 CDS_finder.py --in <FILE> --out <FILE>

Mandatory:
--in     STR     Input file
--out    STR     Output file

Optional:
--len    INT     Minimal CDS length[100]
--atg    STR     Enforces ATG at start [on]
```

`--in` specifies a FASTA file that contains nucleotide sequences.

`--out` specifies a FASTA file that contains the identified coding sequences. Only the longest sequence per input sequence is returned. Sequences below a minimal length threshold can be excluded.

`--len` specifies the minimal CDS length cutoff. Default: 100.

`--atg` specifies whether a CDS must start with an ATG. Possible values are 'on' and 'off'. Default: 'on'.


## References
This repository.
