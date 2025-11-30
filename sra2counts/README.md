# SRA2COUNTS

This README explains how to automatically download FASTQ files from the Sequence Read Archive and to process them with kallisto for the generation of count tables. The script can be started multiple times with the same input and output settings to allow running different SRA data sets in parallel. However, the first job needs some time to prepare the kallisto index which will then be used by all jobs.


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

`--in` specifies a text file with one SRA ID per line.

`--cds` specifies a CDS file that provides the reference sequences for the kallisto analysis.

`--out` specifies the output directory. FASTQ files will be stored here and all temporary files will be stored in subfolders.

`--fastqdump` specifies the fastq-dump path. Default: fastq-dump.

`--kallisto` specifies the kallisto path. Default: kallisto.

`--cpus` specifies the number of CPUs to used by kallisto. Default: 10.

`--min` specifies the minimal percentage of reads that need to be assigned to the top100 transcript for valid RNA-seq samples. Default: 10.

`--max` specifies the maximal percentage of reads that need to be assigned to the top100 transcript for valid RNA-seq samples. Default: 80.

`--mincounts` specifies the minimal number of counts for a sample to be considered as valid RNA-seq data sets. Default: 1000000.



## References
Pucker & Iorizzo, 2022: 

