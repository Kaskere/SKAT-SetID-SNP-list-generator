#!/usr/bin/env Rscript

# A gene/SNPs set list generator (e.g. for SKAT)

# Or Yaacov, or.yaacov@mail.huji.ac.il

# Dependencies: R (3.3+), data.table, tidyverse, Bioconductor.
# Bioconductor packages BSgenome, mygene, TxDb.Hsapiens.UCSC.hg19.knownGene (or hg38).
# If not installed, please install using BiocManager::install("mygene")

# Using the script:
# Change the parameters in the first block of code: Input bim file path, output path and file name, List of gene symbols
# Run the rest of the code.

### change paramerters in this block only ###

bimfile= "/user/plink/myfile.bim" #Plink's bim file (or a exuivalent lisgt of all SNPs in the data)
before= 50000 #number of base pairs to include before the gene
after= 50000 #number of basepairs after

genesymbols= c("NSG1", 
               "CYP11B2") #list of genes

output= "setlist.txt" #output file
library("TxDb.Hsapiens.UCSC.hg19.knownGene") #for hg38 use library("TxDb.Hsapiens.UCSC.hg38.knownGene")
db<- TxDb.Hsapiens.UCSC.hg19.knownGene 

### No changes required after this line, run script ###

library(data.table)
library(BSgenome)
library(tidyverse)
library(mygene)

a<- queryMany(genesymbols, scopes="symbol", fields="entrezgene", species="human")
IDs<- a$entrezgene

#Function to find start position
findStart <- function(geneid) {
  start <- toString(start(ranges(genes(db)[which(genes(db)$gene_id == geneid),])))
  return(start)
}

#Find end position
findEnd <- function(geneid) {
  end <- toString(end(ranges(genes(db)[which(genes(db)$gene_id == geneid),])))
  return(end)
}

#Find chromosome
findChr <- function(geneid) {
    chr <- toString(seqnames(genes(db)[which(genes(db)$gene_id == geneid),]))
    chr <- str_remove(chr, 'chr')
  return(chr)
}

# Make a list of gene positions
chro <- c()
starts <- c()
ends <- c()

for (i in 1:length(genesymbols)) {
    id = IDs[i]
    chro<- append(chro, findChr(id))   
    starts<- append(starts, as.numeric(findStart(id)))
    ends<- append(ends, as.numeric(findEnd(id)))
}

bim <- fread(bimfile)

# extract the SNPs from the bim file:
sets<- character()
snps<- character()
for (i in 1:length(genesymbols)) {
    if (chro[i] <= 22) {
        foundsnps<- bim$V2[which((bim$V1 == chro[i]) & (bim$V4 > (starts[i] - before)) & (bim$V4 < (ends[i] + after)))]
     }
    if (chro[i] == "X") {
        foundsnps<- bim$V2[which(((bim$V1 == 23) | (bim$V1 == 25)) & (bim$V4 > (starts[i] - before)) & (bim$V4 < (ends[i] + after)))]
    }
    if (chro[i] == "Y") {
        foundsnps<- bim$V2[which((bim$V1 == 24) & (bim$V4 > (starts[i] - before)) & (bim$V4 < (ends[i] + after)))]
    }
    for (x in 1:length(foundsnps)) {
    sets<- append(sets, genesymbols[i])
    snps<- append(snps, foundsnps[x]) 
    }     
}

setlist<- data.frame(sets, snps)

#Save output:
write.table(setlist, file=output, sep = "\t", quote = FALSE, col.names = F, row.names = F)
