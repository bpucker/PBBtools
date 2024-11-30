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

`--in` specifies the input file with one SRA accession of a RNA-seq dataset per line.

`--out` specifies the output folder. If this folder does not exist, it will be created. A subfolder will be created to store the count tables belonging to indidvidual SRA accessions.

`--ref` specifies a FASTA file that contains the reference sequences for the kallisto analysis. Usually, these are CDS or transcript sequences.

`--tmp` specifies a temp folder for the storage of downloaded .sra and .fastq files. Subfolders are created for each SRA accession. These folders will be deleted once the process is completed.

`--kallisto` specifies the full path to kallisto. This is only required if a specific kallisto version is to be used or if kallisto is not globally installed on the system. Default: kallisto.

`--fastqdump` specifies  the full path to fastq-dump. This is only required if a specific fastq-dump version is to be used or if fastq-dump is not globally installed on the system. Default: fastq-dump.

`--prefetch` specifies  the full path to prefetch. This is only required if a specific prefetch version is to be used or if prefetch is not globally installed on the system. Default: prefetch.

`--cpus` specifies the number of threads for kallisto. Default: 2.


## Reference
This repository.
