library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
library(ggnewscale)
library(ggridges)
library(europepmc)

# SET THE DESIRED ORGANISM HERE
organism = "org.At.tair.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
organism = org.At.tair.db
# reading in data from deseq2
df = read.csv("all.tsv",sep = "\t", header=TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$GeneID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "TAIR",
             minGSSize = 1, 
             maxGSSize = 1000, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

pwgse <-pairwise_termsim(gse)
emapplot(pwgse, showCategory = 10)

# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)

ridgeplot(gse) + labs(x = "enrichment distribution")

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

terms <- gse$Description[1:3]
pmcplot(terms, 2010:2018, proportion=FALSE)
