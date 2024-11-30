# FASTQ_processor

This script enables the download of RNA-seq data from the Sequence Read Archive (FASTQ files) and processing with kallisto.

## Usage



```
Usage
python3 FASTQ_processor.py --in <FILE> --out <DIR> --ref <FILE> --tmp <DIR>

Mandatory:
--in        STR     Input file with SRA accessions
--out       STR     Output folder
--ref       STR     FASTA file with reference sequences
--tmp       STR     Temp folder

Optional:
--kallisto   STR     Path to kallisto [kallisto]
--fastqdump  STR     Path to fastq-dump [fastq-dump]
--prefetch   STR     Path to prefetch [prefetch]
--cpus       INT     Number of threads to use for kallisto [2]
```

`--in` specifies the input file.

	

## Reference
This repository.
