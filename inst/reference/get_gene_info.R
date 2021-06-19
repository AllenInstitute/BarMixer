library(data.table)

# 10x Genomics GTF file for Reference 3.0.0 (November 19, 2018) for GRCh38
gtf <- fread("C:/datasets/genes.gtf",
             skip = 5)
names(gtf) <- c("chr","source","type","start","end","V6","strand","V8","meta")
gtf <- gtf[type == "gene",]

gtf$id <- sub("\";.+","",
              sub("gene_id \"","",gtf$meta))
gtf$name <- sub("\";.+","",
                sub(".+gene_name \"","",gtf$meta))
gtf$biotype <- sub("\"","",
                   sub(".+gene_biotype \"","",gtf$meta))

gtf$chr <- paste0("chr",gtf$chr)

gtf$meta <- NULL
gtf$source <- NULL
gtf$type <- NULL
gtf$V6 <- NULL
gtf$V8 <- NULL

gene_meta <- setcolorder(gtf,
                         c('chr', 'start', 'end', 'strand', 'id', 'name', 'biotype'))

fwrite(gene_meta,
       "inst/reference/GRCh38_10x_gene_metadata.csv.gz")

chrM_genes <- gene_meta[chr == "chrMT",]

fwrite(chrM_genes,
       "inst/reference/GRCh38_10x_chrM_gene_metadata.csv.gz")
