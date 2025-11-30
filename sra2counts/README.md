# SRA2COUNTS

This README explains how to automatically download FASTQ files from the Sequence Read Archive and to process them with kallisto for the generation of count tables.


## Usage ##

```
Usage
python3 reads2counts.py --in <FILE> --cds <FILE> --out <OUTPUT_DIR>

Mandatory:
--in     STR     Input file containing SRA IDs
--cds    STR     CDS file for reference
--out    STR     Output directory

Optional:
--fastqdump   STR   fastq-dump path [fastq-dump]
--kallisto <FULL_PATH_TO_KALLISTO>[kallisto]
--cpus <NUMBER_OF_CPUS_TO_USE>[10]
--min <MIN_PERCENT_EXPRESSION_ON_TOP100>[10]
--max <MAX_PERCENT_EXPRESSION_ON_TOP100>[80]
--mincounts <MIN_READ_NUMBER>[1000000]
```

`--in` specifies a FASTA file/directory that contains coding sequences.



## References
Pucker & Iorizzo, 2022: 

