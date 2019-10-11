# SKAT-SetID-SNP-list-generator
Takes a Plinks BIM file and a list of gene symbols, and generates a SKAT SetID list file

#### Dependencies: 
* R (3.3+), data.table, tidyverse, Bioconductor.
* Bioconductor packages BSgenome, mygene, TxDb.Hsapiens.UCSC.hg19.knownGene (or hg38).
* If not installed, please install using BiocManager::install("mygene")

#### Using the script:
Change the parameters in the first block of code:
* Input bim file path, output path and file name
* List of gene symbols
* Run the rest of the code.

or.yaacov@mail.huji.ac.il
