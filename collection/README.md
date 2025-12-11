# Documentation of script usage

Please see the usage in the individual scripts for details.

## construct_anno.py
This script annotates novel polypeptide sequences by comparison against a reference (typically _A. thaliana_). Annotation text from the reference sequences is transferred to the novel sequences. Therefore, reciprocal best BLAST hits (RBHs) are identified between the novel sequences and the reference sequences. These are likely to be orthologs. If no orthologs are detected, the best hit for each novel sequence is utilized.
Temporary files include the files beloning to BLAST databases that are created for the input FASTA file and the reference FASTA file to enable the RBH analysis. There are also temporary files with the BLAST hits for ortholog detection via RBHs. RBHs are written into a dedicated text file.


```
Usage:
  python3 construct_anno.py

Mandatory:
  --in      STR     FASTA file containing peptide sequences for annotation.
  --out     STR     Output folder
  --ref     STR     Reference sequence FASTA file
  --anno    STR     Annotation file (TSV)

```


`--in` specifies the input FASTA file.

`--out` specifies folder. Will be created if it does not exist already.

`--ref` specifies the reference sequence FASTA file.

`--anno` specifies the annotation file (TSV).
