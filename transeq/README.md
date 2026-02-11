# Translation of CDS into peptide sequence
This script reads all sequences from a given FASTA file and converts them into the corresponding peptide sequences. It translates unknown/ambiguous codons (e.g., containing N) to 'X'. If --internal-stop-to-x is provided, internal stop codons (represented by '__\*__') are converted to 'X' (a final '*' is kept).

## Usage ##

```
Usage
python3 transeq.py --in <FILE|INPUT_DIR> --out <FILE|OUTPUT_DIR>

Mandatory:
--in     STR     Input file or directory containing (multiple) FASTA files
--out    STR     Output file or directory

Optional:
--internal-stop-to-x  Preset  If provided, internal stop codons converted to 'X'
```

`--in` specifies a FASTA file/directory that contains coding sequences.

`--out` specifies a FASTA file/directory that contains the generated peptide sequences.

