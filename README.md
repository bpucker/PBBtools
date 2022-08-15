# PBBTools
collection of mini tools for the BioinfToolServer

# Identification of Open Reading Frame
This script screens given nucleotide sequences for the longes Open Reading Frame (ORF) / coding sequence (CDS) that could encode a peptide sequence. The longest detected ORF/CDS per input sequence is reported.

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


# Translation of CDS into peptide sequence
This script reads all sequences from a given FASTA file and converts them into the corresponding peptide sequences.

## Usage ##

```
Usage
python3 transeq.py --in <FILE> --out <FILE>

Mandatory:
--in     STR     Input file
--out    STR     Output file
```

`--in` specifies a FASTA file that contains coding sequences.

`--out` specifies a FASTA file that contains the generated peptide sequences.
