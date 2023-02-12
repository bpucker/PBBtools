# Dotplotter

Visualizing similarity of two sequences through a dot plot.


## Usage ##

```
Usage:
python3 dotplotter.py --in1 <FILE> --in2 <FILE> --out <FILE>
or
python3 dotplotter.py --seq1 <STR> --seq2 <STR> --out <FILE>

Option1:
--in1     STR     Input FASTA file1
--in2     STR     Input FASTA file2
--out     STR     Figure output file

Option2:
--seq1    STR     Sequence1
--seq2    STR     Sequence2
--out     STR     Figure output file

Optional:
--kmer    INT     k-mer length[31]
--name1   STR     Sequence1 name
--name2   STR     Sequence2 name
```

`--in1` specifies a FASTA file that contains the first nucleotide sequence. The sequence name from this file can be used as label if no other name is specified via `--name1`.

`--in2` specifies a FASTA file that contains the second nucleotide sequence. The sequence name from this file can be used as label if no other name is specified via `--name2`.

`--seq1` specifies the first input sequence. The characters will be turned into upper case characters to avoid mismatches due to case sensitivity.

`--seq2` specifies the second input sequence. The characters will be turned into upper case characters to avoid mismatches due to case sensitivity.

`--out` specifies the figure output file. The file extension determines the type of figure file that is created. Only types supported by Python and the operating system are possible. Usually PNG, JPG, SVG, and PDF are supported.

`--kmer` specifies the size of k-mers to use for identification of similarities between both sequences. Smaller k-mer sizes will lead to more matches, but larger values will make the comparison more specific. Default: 31.

`--name1` specifies the name of sequence 1 that will be displayed as axis label in the dotplot. This argument overwrites the sequence name that might have been retrieved from sequence file 1 if the sequence was provided as FASTA file.

`--name2` specifies the name of sequence 2 that will be displayed as axis label in the dotplot. This argument overwrites the sequence name that might have been retrieved from sequence file 2 if the sequence was provided as FASTA file.


## References

This repository.
