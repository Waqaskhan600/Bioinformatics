#!/usr/bin/env Rscript
# differential_gene_expression_analysis.R
require(optparse)

#Parse arguments from command line
options <- list(
 make_option(c(" - input_count_data_path"), action = "store", type = "character", help="Path to the raw count data."),
 make_option(c(" - metadata_path"), action = "store", type = "character", help="Path to the metadata."),
 make_option(c(" - outdir"), action = "store", default = "data", type = "character", help="Output directory for storing differential gene expression analysis results."),
 make_option(c(" - factor_level1"), action = "store", default = "Type", help="Class level 1 column name."),
 make_option(c(" - factor_level2"), action = "store", default = "Type", help="Class level 2 column name.")
)
arguments <- parse_args(OptionParser(option_list = options))

# Check if mandatory arguments are provided
if (is.null(arguments$input_count_data_path) || is.null(arguments$metadata_path)) {
 stop("Error: Mandatory arguments - input_count_data_path and - metadata_path are required.")
}

# create variable handles for parsed parameters
input_count_data_path <- arguments$input_count_data_path
metadata_path <- arguments$metadata_path
outdir <- arguments$outdir
factor_level1 <- arguments$factor_level1
factor_level2 <- arguments$factor_level2

# Load DESeq2 library
library(DESeq2);
library(dplyr);
library(pheatmap)
library(ggplot2)
library("RColorBrewer")
library(ggrepel)

# create output directory
outdir <- paste0(outdir, "/dge")
dir.create(outdir)

# Load count data
count_data <- read.table(input_count_data_path, header = TRUE, row.names = 1, sep = ",")
colnames(count_data)

# load metadata
metadata <- read.table(metadata_path, header = TRUE, row.names = 1, sep = ",")
all(rownames(metadata) %in% colnames(count_data))
all(rownames(metadata) == colnames(count_data))

# set factor levels
factors <- factor(metadata[[factor_level1]])
metadata[[factor_level1]] <- factors

# Create a deseq object and import the count data and metatdata
dds <- DESeqDataSetFromMatrix(countData = round(count_data),
 colData = metadata,
 design = ~ Type)

# set reference for group factors
dds[[factor_level1]] <- relevel(dds[[factor_level1]], ref = "healthy control")

# filter out low count genes
# keep genes with at least N counts >= 10, where N = size of smallest group
keep <- rowSums(counts(dds) >=10) >= min(table(metadata[[factor_level1]]))
dds <- dds[keep,]

# Perform Statistical Test using Deseq
dds <- DESeq(dds, test = "Wald")

# get results of Deseq
deseq_result <- results(dds, cooksCutoff = FALSE, independentFiltering = FALSE)

# process results
deseq_result_df <- as.data.frame(deseq_result)
deseq_result_df <- cbind(GeneName = row.names(deseq_result_df), deseq_result_df)

# save results
write.table(deseq_result_df, file = paste0(outdir, "/Deseq2-results-all.tsv"), row.names = F, sep="\t")

# extract genes with padj < 0.05 and log2Foldchange <= -1 or >= 1
deg <- subset(deseq_result, padj<0.05 & abs(log2FoldChange) >=1)

# sort by adjusted p-values
deg <- deg[order(deg$padj),]

# save filtered results
write.table(deg, file = paste0(outdir, "/Deseq2-results-filtered.tsv"), row.names = F, sep = "\t")

#########################################
# Gene expression data visualization
#########################################
################################################################################
# Start MA Plot
resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm") # remove noise
# Save the plot as a PDF/PNG file
pdf(paste0(outdir, "/MA_plot.pdf"), width = 8, height = 6)
plotMA(deseq_result, ylim = c(-2, 2))
dev.off()
png(paste0(outdir, "/MA_plot.png"), width = 480, height = 360)
plotMA(deseq_result, ylim = c(-2, 2))
dev.off()
pdf(paste0(outdir, "/MA_plot_denoised.pdf"), width = 8, height = 6)
plotMA(resLFC, ylim=c(-2,2))
dev.off()
png(paste0(outdir, "/MA_plot_denoised.png"), width = 480, height = 360)
plotMA(resLFC, ylim=c(-2,2))
dev.off()
# End

################################################################################
################################################################################
# Start Dispersion plot
pdf(paste0(outdir, "/Dispersion_plot.pdf"), width = 8, height = 6)
plotDispEsts(dds)
dev.off()
png(paste0(outdir, "/Dispersion_plot.png"), width = 480, height = 360)
plotDispEsts(dds)
dev.off()
# End

################################################################################
################################################################################
# Start volcano plots # set groups
old.pal <- palette(c("#FF0000", "#0000FF", "#808080")) # set colors
par(mar=c(4,4,2,1), cex.main=1.5) # set margin size
pdf(paste0(outdir, "/Volcano_plot.pdf"), width = 8, height = 6)

#plot values
plot(deseq_result$log2FoldChange, -log10(deseq_result$padj), col = "darkgray",
 xlab="log2FoldChange", ylab="-log10(Padj)", pch=20, cex=0.5)
with( subset(deseq_result, padj <0.05 & abs(log2FoldChange) >=1),
 points(log2FoldChange, -log10(padj), pch =20, col=(sign(log2FoldChange) +3)/2, cex =1))
legend("bottomleft", title=paste("Padj<", 0.05,sep=""),
 legend=c("down regulated", "up regulated", "unchanged"), pch=20, col=1:3)
dev.off()
resLFC_df <- as.data.frame(deseq_result_df)
resLFC_df$diffexpressed <- "no change"
resLFC_df$diffexpressed[resLFC_df$log2FoldChange>1 & resLFC_df$padj<0.05] <- "up regulated"
resLFC_df$diffexpressed[resLFC_df$log2FoldChange< (-1) & resLFC_df$padj<0.05] <- "down regulated"
resLFC_df$diffexpressed[resLFC_df$padj>=0.05] <- "insignificant"
resLFC_df$delabel <- NA
volcano2 <- ggplot(data = resLFC_df, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed,label=GeneName))+
 geom_point()+
 theme_minimal()+
 geom_text_repel()+
 scale_color_manual(values=c("red","black","darkgray","blue"))+
 theme(text=element_text(size=10))
ggsave(filename = paste0(outdir, "/Volcano_plot2.png"), plot = volcano2 , bg = "white",device = "png", width = 8, height = 6)
ggsave(filename = paste0(outdir, "/Volcano_plot2.pdf"), plot = volcano2 , device = "pdf", width = 8, height = 6)
# End
#############################################################################