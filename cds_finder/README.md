## Usage ##

```
Usage
python3 CDS_finder.py --in <FILE> --out <FILE>

Mandatory:
--in     STR     Input file
--out    STR     Output file

Optional:
--len    INT     Minimal CDS length[100]
```

`--in` specifies a FASTA file that contains nucleotide sequences.

`--out` specifies a FASTA file that contains the identified coding sequences. Only the longest sequence per input sequence is returned. Sequences below a minimal length threshold can be excluded.

`--len` specifies the minimal CDS length cutoff. Default: 100.
