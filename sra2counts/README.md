This README explains how to use the tool:

python3 ./SRA2counts.py \
--sra ./SRA_samples.txt \
--cds ./Araport11_genes.201606.cds.repr_MOD.fasta \
--out ./test/ \
--downloader./automatic_SRA_download.py \
--pipeline ./kallisto_pipeline3.py \
--merger ./merge_kallisto_output3.py \
--cleaner ./filter_RNAseq_samples.py
