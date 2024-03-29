### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###

### DGE analysis for Arabidopsis RNA-Seq samples in Applied Genome Research course ###

library("DESeq2")

# --- loading sampleTable --- #
csvfile <- "./clean_data/clean_sample_table.txt"
sampleTable <- read.csv(csvfile, row.names=1, sep="\t")
sampleTable$genotype <- as.factor( sampleTable$genotype )

summary(sampleTable)

# --- loading the data matrix --- #
count_data_file <- "./clean_data/clean_data_matrix.txt"
countdata <- read.csv(count_data_file,row.names=1, header=T, sep="\t")
summary(countdata)

# --- construction of DESeqDataSet --- #
ddsMat <- DESeqDataSetFromMatrix( countData=countdata, colData=sampleTable, design= ~ genotype )
nrow(ddsMat)

# -- removal of not or low expressed genes --- #
dds <- ddsMat[ rowSums(counts(ddsMat)) > 1, ]
nrow(dds)

# --- plot PCA in R studio --- #
rld <- rlog(dds)
ramp <- 1:2/2
cols <- c( rgb(ramp, 0, 0), rgb(0, ramp, 0), rgb(ramp, 0, ramp), rgb(ramp, 0, ramp) )
print ( plotPCA( rld, intgroup=c( "genotype" )  ) )

# --- differential expression analysis --- #
dds <- DESeq(dds)
res <- results(dds)
summary(res)

# --- investigate differentially expressed genes --- #
res.05 <- results( dds, alpha=.05  )
table(res.05$padj < .05  )

write.table( res.05, "./20201029_results_table.txt" )

# --- print session information --- #
sessionInfo()

