# SRA2COUNTS

This README explains how to automatically download FASTQ files from the Sequence Read Archive and to process them with kallisto for the generation of count tables.



					python3 reads2counts.py (""" +  __version__	+ """)
					--in <FILE_WITH_ONE_RUN_ID_PER_LINE>
					--cds <CDS_FILE>
					--out <OUTPUT_DIRECTORY>
					
					optional:
					--fastqdump <FULL_PATH_TO_FASTQ_DUMP>[fastq-dump]
					--kallisto <FULL_PATH_TO_KALLISTO>[kallisto]
					--cpus <NUMBER_OF_CPUS_TO_USE>[10]
					--min <MIN_PERCENT_EXPRESSION_ON_TOP100>[10]
					--max <MAX_PERCENT_EXPRESSION_ON_TOP100>[80]
					--mincounts <MIN_READ_NUMBER>[1000000]
