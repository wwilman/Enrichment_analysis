#BiocManager::install("clusterProfiler", version = "3.8")
#BiocManager::install("pathview")
#install.packages("wordcloud")
library(clusterProfiler)
library(wordcloud)

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

# Exctract significant results (padj < 0.05)
sig_genes_df = subset(df, padj < 0.05)

# From significant results, we want to filter on log2fold change
genes <- sig_genes_df$log2FoldChange

# Name the vector
names(genes) <- sig_genes_df$GeneID

# omit NA values
genes <- na.omit(genes)

# filter on min log2fold change (log2FoldChange > 2)
genes <- names(genes)[abs(genes) > 2]

go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'TAIR',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
library(enrichplot)
library(ggupset)
upsetplot(go_enrich)

wcdf<-read.table(text=go_enrich$GeneRatio, sep = "/")[1]
wcdf$term<-go_enrich[,2]
wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(1, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)

dotplot(go_enrich)

pwgse <-pairwise_termsim(go_enrich)
emapplot(pwgse)

goplot(go_enrich, showCategory = 10)

# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(go_enrich, categorySize="pvalue", foldChange=gene_list)
