[![DOI](https://zenodo.org/badge/524938056.svg)](https://doi.org/10.5281/zenodo.7720916)

# PBBtools: bioinformatic tools developed by PuckerLab
This is a collection of tools and useful scripts developed by members of the [PuckerLab](https://www.pbb.uni-bonn.de). Many of these tools are freely available through our [BioinfToolServer](https://www.pbb-tools.de). Please [get in touch](https://www.izmb.uni-bonn.de/en/molecular-plant-sciences/contact) if you have any questions or suggestions for additional features.



# KIPEs: Knowledge-based Identification of Pathway Enzymes
**Description**: This tool facilitates the identification of genes/proteins involved in a particular pathway. KIPEs was initially developed for the annotation of the flavonoid biosynthesis, but can also be applied to other pathways.

**Link**: [KIPEs repository](https://github.com/bpucker/KIPEs) and [KIPEs webserver](https://pbb-tools.de/KIPEs/)

**References**: 
Rempel A., Choudhary N. and Pucker B. (2023). KIPEs3: Automatic annotation of biosynthesis pathways. PLOS ONE 18(11): e0294342. doi: [10.1371/journal.pone.0294342](https://doi.org/10.1371/journal.pone.0294342).

Pucker, B.; Reiher, F.; Schilbert, H.M. Automatic Identification of Players in the Flavonoid Biosynthesis with Application on the Biomedicinal Plant Croton tiglium. Plants 2020, 9, 1103. doi:[10.3390/plants9091103](https://doi.org/10.3390/plants9091103).


# MYB_annotator: Automatic identification and annotation of MYB transcription factors
**Description**: This tool facilitates the identification of MYB transcription factors in a given species. Identified MYBs are assigned to well characterized orthologs and a functional annotation is transferred.

**Link**: [MYB_annotator repository](https://github.com/bpucker/MYB_annotator) and [MYB_annotator webserver](https://pbb-tools.de/MYB_annotator/)

**Reference**: Pucker, B. Automatic identification and annotation of MYB gene family members in plants. BMC Genomics 23, 220 (2022). doi:[10.1186/s12864-022-08452-5](https://doi.org/10.1186/s12864-022-08452-5).



# bHLH_annotator: Automatic identification and annotation of bHLH transcription factors
**Description**: This tool facilitates the identification of bHLH transcription factors in a given species. Identified bHLHs are assigned to well characterized orthologs and a functional annotation is transferred.

**Link**: [bHLH_annotator repository](https://github.com/bpucker/bHLH_annotator) and [bHLH_annotator webserver](https://pbb-tools.de/bHLH_annotator/)

**Reference**: Thoben C. and Pucker B. (2023). Automatic annotation of the bHLH gene family in plants. BMC Genomics 24, 780 (2023). doi: [10.1186/s12864-023-09877-2](https://doi.org/10.1186/s12864-023-09877-2).





# CoExp: Co-expression analysis
**Description**: This tool performs a pairwise co-expression analysis for a given set of genes of interest. All co-expressed genes are identified for this set.

**Link**: [CoExp repository](https://github.com/bpucker/CoExp) and [CoExp webserver](https://pbb-tools.de/CoExp/)

**References**: 
Pucker B, Iorizzo M (2023) Apiaceae FNS I originated from F3H through tandem gene duplication. PLOS ONE 18(1): e0280155. doi:[10.1371/journal.pone.0280155](https://doi.org/10.1371/journal.pone.0280155).

Pucker B., Walker-Hale N., Yim W.C., Cushman J.C., Crumm A., Yang Y., Brockington S. (2024). Evolutionary blocks to anthocyanin accumulation and the loss of an anthocyanin carrier protein in betalain-pigmented Caryophyllales. New Phytol, 241: 471-489. doi:[10.1111/nph.19341](https://doi.org/10.1111/nph.19341).



# HeatmapPlotter: Gene expression visualization
**Description**: This script generates a heatmap for a selection of genes and RNA-seq samples.

**Link**: [HeatmapPlotter repository](https://github.com/bpucker/PBBtools/blob/main/heatmapper/README.md)

**References**: This repository.


# CDS_finder: Identification of Open Reading Frame
**Description**: This script screens given nucleotide sequences for the longes Open Reading Frame (ORF) / coding sequence (CDS) that could encode a peptide sequence. The longest detected ORF/CDS per input sequence is reported.

**Link**: [CDS_finder repository](https://github.com/bpucker/PBBTools/blob/main/cds_finder/README.md)

**Reference**: This repository.



# TranSeq: Translation of CDS into peptide sequence
**Description**: This script reads all sequences from a given FASTA file and converts them into the corresponding peptide sequences.

**Link**: [transeq repository](https://github.com/bpucker/PBBTools/blob/main/transeq/README.md)

**Reference**: This repository.



# SeqEx: Extracting sequence of interest from FASTA file
**Description**: This script extracts a specified region on a given sequence from a FASTA file.

**Link**: [seqex repository](https://github.com/bpucker/PBBtools/blob/main/seqex/README.md)

**Reference**: Pucker B, Schilbert H, Schumacher SF. Integrating Molecular Biology and Bioinformatics Education. Journal of Integrative Bioinformatics. 2019;16(3): 20190005. doi:[10.1515/jib-2019-0005](https://doi.org/10.1515/jib-2019-0005).



# Dotplotter: Visualization of sequence similarity in dot plot
**Description**: This script compares to sequences and visualizes matches by dots in a 2D plot.

**Link**: [dotplotter repository](https://github.com/bpucker/PBBtools/blob/main/dotplotter/README.md)

**Reference**: This repository.



# MGSE: Mapping-based Genome Size Estimation
**Description**: This script estimates the size of a plant genome based on the coverage in a read mapping.

**Link**: [MGSE repository](https://github.com/bpucker/MGSE)

**Reference**: Natarajan, S., Gehrke, J. & Pucker, B. Mapping-based genome size estimation. BMC Genomics 26, 482 (2025). doi: [10.1186/s12864-025-11640-8](https://doi.org/10.1186/s12864-025-11640-8).



# NAVIP: Neighborhood-Aware Variant Impact Predictor
**Description**: This script predicts the functional impact of sequence variants on the surrounding gene. NAVIP considers all variants in a gene simultaneously to account for their interactions.

**Link**: [NAVIP repository](https://github.com/bpucker/NAVIP)

**Reference**: Baasner J-S, Rempel A, Howard D, Pucker B (2025) NAVIP: Unraveling the influence of neighboring small sequence variants on functional impact prediction. PLoS Comput Biol 21(2): e1012732. doi: [10.1371/journal.pcbi.1012732](https://doi.org/10.1371/journal.pcbi.1012732).



# Long Read Walker
**Description**: This script closes the gaps between two contigs in a long read assembly by the extension of overlaps.

**Link**: [LongReadWalker repository](https://github.com/bpucker/LongReadWalker)

**Reference**: Sielemann, K., Pucker, B., Orsini, E. et al. Genomic characterization of a nematode tolerance locus in sugar beet. BMC Genomics 24, 748 (2023). doi: [10.1186/s12864-023-09823-2](https://doi.org/10.1186/s12864-023-09823-2).


# Documentation of methods
[Dry lab methods](https://github.com/bpucker/pbb-drylab)

[Wet lab methods](https://github.com/bpucker/pbb-wetlab)



