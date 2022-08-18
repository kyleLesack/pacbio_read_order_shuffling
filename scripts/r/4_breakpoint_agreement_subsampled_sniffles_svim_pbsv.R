library(ggplot2)
library(tidyverse)
library(patchwork)
setwd("/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling")

pbsv_pbmm2_10X.df <- read.delim(file = '4_results/subsampled/10X/breakpoint_agreement/pbsv/ggplot/pbsv-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
pbsv_pbmm2_20X.df <- read.delim(file = '4_results/subsampled/20X/breakpoint_agreement/pbsv/ggplot/pbsv-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
pbsv_pbmm2_40X.df <- read.delim(file = '4_results/subsampled/40X/breakpoint_agreement/pbsv/ggplot/pbsv-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
pbsv_pbmm2_60X.df <- read.delim(file = '4_results/subsampled/60X/breakpoint_agreement/pbsv/ggplot/pbsv-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)

pbsv_discordance.df <- rbind(pbsv_pbmm2_10X.df, pbsv_pbmm2_20X.df, pbsv_pbmm2_40X.df, pbsv_pbmm2_60X.df)
pbsv_discordance.df
col_order <- c("Depth", "Aligner", "Agreement",  "Count")
pbsv_discordance_reordered.df <- pbsv_discordance.df[, col_order]
pbsv_discordance_reordered.df <- pbsv_discordance_reordered.df[- grep("breakpoints", pbsv_discordance_reordered.df$Agreement),]
pbsv_discordance_reordered.df

#ggplot(pbsv_discordance_reordered.df, aes(x=Depth, y = Count, fill=Aligner)) + geom_col( position=position_dodge())
p_pbsv <- ggplot(pbsv_discordance_reordered.df, aes(x=Depth, y = Count)) + geom_col( position=position_dodge()) + facet_grid(~ Aligner, scales = "free_x", space = "free_x", switch = "x") + ylim(0, 100)
p_pbsv <- p_pbsv + xlab("pbsv") + ylab("Discordant Breakpoints (%)") 
#p_pbsv <- p_pbsv + xlab("pbsv") + theme(axis.title.y = element_blank())
#p_pbsv <- p_pbsv + xlab("pbsv") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())
#p_pbsv <- p_pbsv + xlab("Sequencing Depth") + ylab("Discordant Breakpoints (%)") + ggtitle("Impact of Sequencing Depth on pbsv Discordant Breakpoints") + theme(plot.title = element_text(hjust = 0.5)) #+ theme_bw()

sniffles_minimap2_10X.df <- read.delim(file = '4_results/subsampled/10X/breakpoint_agreement/sniffles/ggplot/sniffles-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_ngmlr_10X.df <- read.delim(file = '4_results/subsampled/10X/breakpoint_agreement/sniffles/ggplot/sniffles-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_pbmm2_10X.df <- read.delim(file = '4_results/subsampled/10X/breakpoint_agreement/sniffles/ggplot/sniffles-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_minimap2_20X.df <- read.delim(file = '4_results/subsampled/20X/breakpoint_agreement/sniffles/ggplot/sniffles-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_ngmlr_20X.df <- read.delim(file = '4_results/subsampled/20X/breakpoint_agreement/sniffles/ggplot/sniffles-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_pbmm2_20X.df <- read.delim(file = '4_results/subsampled/20X/breakpoint_agreement/sniffles/ggplot/sniffles-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_minimap2_40X.df <- read.delim(file = '4_results/subsampled/40X/breakpoint_agreement/sniffles/ggplot/sniffles-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_ngmlr_40X.df <- read.delim(file = '4_results/subsampled/40X/breakpoint_agreement/sniffles/ggplot/sniffles-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_pbmm2_40X.df <- read.delim(file = '4_results/subsampled/40X/breakpoint_agreement/sniffles/ggplot/sniffles-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_minimap2_60X.df <- read.delim(file = '4_results/subsampled/60X/breakpoint_agreement/sniffles/ggplot/sniffles-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_ngmlr_60X.df <- read.delim(file = '4_results/subsampled/60X/breakpoint_agreement/sniffles/ggplot/sniffles-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_pbmm2_60X.df <- read.delim(file = '4_results/subsampled/60X/breakpoint_agreement/sniffles/ggplot/sniffles-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)

sniffles_discordance.df <- rbind(sniffles_minimap2_10X.df, sniffles_minimap2_20X.df, sniffles_minimap2_40X.df, sniffles_minimap2_60X.df,sniffles_ngmlr_10X.df, sniffles_ngmlr_20X.df, sniffles_ngmlr_40X.df, sniffles_ngmlr_60X.df,sniffles_pbmm2_10X.df, sniffles_pbmm2_20X.df, sniffles_pbmm2_40X.df, sniffles_pbmm2_60X.df)
sniffles_discordance.df
col_order <- c("Depth", "Aligner", "Agreement",  "Count")
sniffles_discordance_reordered.df <- sniffles_discordance.df[, col_order]
sniffles_discordance_reordered.df <- sniffles_discordance_reordered.df[- grep("breakpoints", sniffles_discordance_reordered.df$Agreement),]
sniffles_discordance_reordered.df

#ggplot(sniffles_discordance_reordered.df, aes(x=Depth, y = Count, fill=Aligner)) + geom_col( position=position_dodge())
p_sniffles <- ggplot(sniffles_discordance_reordered.df, aes(x=Depth, y = Count)) + geom_col( position=position_dodge()) + facet_grid(~ Aligner, scales = "free_x", space = "free_x", switch = "x") + ylim(0, 100)
p_sniffles <- p_sniffles + xlab("Sniffles") + ylab("Discordant Breakpoints (%)") 

svim_minimap2_10X.df <- read.delim(file = '4_results/subsampled/10X/breakpoint_agreement/svim/qual_15/ggplot/svim-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_ngmlr_10X.df <- read.delim(file = '4_results/subsampled/10X/breakpoint_agreement/svim/qual_15/ggplot/svim-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_pbmm2_10X.df <- read.delim(file = '4_results/subsampled/10X/breakpoint_agreement/svim/qual_15/ggplot/svim-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_minimap2_20X.df <- read.delim(file = '4_results/subsampled/20X/breakpoint_agreement/svim/qual_15/ggplot/svim-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_ngmlr_20X.df <- read.delim(file = '4_results/subsampled/20X/breakpoint_agreement/svim/qual_15/ggplot/svim-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_pbmm2_20X.df <- read.delim(file = '4_results/subsampled/20X/breakpoint_agreement/svim/qual_15/ggplot/svim-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_minimap2_40X.df <- read.delim(file = '4_results/subsampled/40X/breakpoint_agreement/svim/qual_15/ggplot/svim-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_ngmlr_40X.df <- read.delim(file = '4_results/subsampled/40X/breakpoint_agreement/svim/qual_15/ggplot/svim-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_pbmm2_40X.df <- read.delim(file = '4_results/subsampled/40X/breakpoint_agreement/svim/qual_15/ggplot/svim-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_minimap2_60X.df <- read.delim(file = '4_results/subsampled/60X/breakpoint_agreement/svim/qual_15/ggplot/svim-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_ngmlr_60X.df <- read.delim(file = '4_results/subsampled/60X/breakpoint_agreement/svim/qual_15/ggplot/svim-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_pbmm2_60X.df <- read.delim(file = '4_results/subsampled/60X/breakpoint_agreement/svim/qual_15/ggplot/svim-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)

svim_discordance.df <- rbind(svim_minimap2_10X.df, svim_minimap2_20X.df, svim_minimap2_40X.df, svim_minimap2_60X.df,svim_ngmlr_10X.df, svim_ngmlr_20X.df, svim_ngmlr_40X.df, svim_ngmlr_60X.df,svim_pbmm2_10X.df, svim_pbmm2_20X.df, svim_pbmm2_40X.df, svim_pbmm2_60X.df)
col_order <- c("Depth", "Aligner", "Agreement",  "Count")
svim_discordance_reordered.df <- svim_discordance.df[, col_order]
svim_discordance_reordered.df <- svim_discordance_reordered.df[- grep("breakpoints", svim_discordance_reordered.df$Agreement),]

#ggplot(svim_discordance_reordered.df, aes(x=Depth, y = Count, fill=Aligner)) + geom_col( position=position_dodge())
p_svim <- ggplot(svim_discordance_reordered.df, aes(x=Depth, y = Count)) + geom_col( position=position_dodge()) + facet_grid(~ Aligner, scales = "free_x", space = "free_x", switch = "x") + ylim(0, 100)
p_svim <- p_svim + xlab("SVIM")  + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())

combined_plots <- p_sniffles + p_svim + p_pbsv + plot_layout(nrow = 2, byrow = TRUE)
combined_plots <- combined_plots + plot_annotation(title = 'Impact of Sequencing Depth on Discordant Breakpoints', theme = theme(plot.title = element_text(hjust = 0.5)))


#ylab <- p_sniffles$labels$y
#p_sniffles$labels$y <- p_svim$labels$y <- " "
#combined_plots

#grid::grid.draw(grid::textGrob(ylab, x = 0.02, rot = 90))

dir.create(file.path("5_plots","subsampled"), recursive = TRUE)
ggsave("4_breakpoint_agreement_sniffles_svim_pbsv_subsampled.png", plot = combined_plots, device = "png", path = "5_plots/subsampled", scale = 1, dpi = 300, limitsize = TRUE, width = 2605, height = 1715, units = "px")

# Plot with counts on top of bars

p_pbsv_counts <- ggplot(pbsv_discordance_reordered.df, aes(x=Depth, y = Count, group = Depth)) + geom_col( position=position_dodge()) + facet_grid(~ Aligner, scales = "free_x", space = "free_x", switch = "x") + ylim(0, 100) + geom_text(aes(x=Depth, y = Count, label = Count, group = Depth), position = position_dodge(width = 1),vjust = -0.5, size = 2.5)
p_pbsv_counts <- p_pbsv_counts + xlab("pbsv") + ylab("Discordant Breakpoints (%)") 

p_sniffles_counts <- ggplot(sniffles_discordance_reordered.df, aes(x=Depth, y = Count, group = Depth)) + geom_col( position=position_dodge()) + facet_grid(~ Aligner, scales = "free_x", space = "free_x", switch = "x") + ylim(0, 100) + geom_text(aes(x=Depth, y = Count, label = Count, group = Depth), position = position_dodge(width = 1),vjust = -0.5, size = 2.5)
p_sniffles_counts <- p_sniffles_counts + xlab("Sniffles") + ylab("Discordant Breakpoints (%)") 

p_svim_counts <- ggplot(svim_discordance_reordered.df, aes(x=Depth, y = Count, group = Depth)) + geom_col( position=position_dodge()) + facet_grid(~ Aligner, scales = "free_x", space = "free_x", switch = "x") + ylim(0, 100) + geom_text(aes(x=Depth, y = Count, label = Count, group = Depth), position = position_dodge(width = 1),vjust = -0.5, size = 2.5)
p_svim_counts <- p_svim_counts + xlab("SVIM")  + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())

combined_plots_counts <- p_sniffles_counts + p_svim_counts + p_pbsv_counts + plot_layout(nrow = 2, byrow = TRUE)
combined_plots_counts <- combined_plots_counts + plot_annotation(title = 'Impact of Sequencing Depth on Discordant Breakpoints', theme = theme(plot.title = element_text(hjust = 0.5)))
ggsave("4_breakpoint_agreement_sniffles_svim_pbsv_subsampled-counts.png", plot = combined_plots_counts, device = "png", path = "5_plots/subsampled", scale = 1, dpi = 300, limitsize = TRUE, width = 2605, height = 1715, units = "px")
