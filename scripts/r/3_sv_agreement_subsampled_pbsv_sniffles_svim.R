library(ggplot2)
library(tidyverse)
library(patchwork)
setwd("/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling")

sniffles_minimap2_10X.df <- read.delim(file = '4_results/subsampled/10X/sv_intersection_agreement/sniffles/ggplot/sniffles-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_ngmlr_10X.df <- read.delim(file = '4_results/subsampled/10X/sv_intersection_agreement/sniffles/ggplot/sniffles-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_pbmm2_10X.df <- read.delim(file = '4_results/subsampled/10X/sv_intersection_agreement/sniffles/ggplot/sniffles-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_minimap2_20X.df <- read.delim(file = '4_results/subsampled/20X/sv_intersection_agreement/sniffles/ggplot/sniffles-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_ngmlr_20X.df <- read.delim(file = '4_results/subsampled/20X/sv_intersection_agreement/sniffles/ggplot/sniffles-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_pbmm2_20X.df <- read.delim(file = '4_results/subsampled/20X/sv_intersection_agreement/sniffles/ggplot/sniffles-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_minimap2_40X.df <- read.delim(file = '4_results/subsampled/40X/sv_intersection_agreement/sniffles/ggplot/sniffles-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_ngmlr_40X.df <- read.delim(file = '4_results/subsampled/40X/sv_intersection_agreement/sniffles/ggplot/sniffles-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_pbmm2_40X.df <- read.delim(file = '4_results/subsampled/40X/sv_intersection_agreement/sniffles/ggplot/sniffles-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_minimap2_60X.df <- read.delim(file = '4_results/subsampled/60X/sv_intersection_agreement/sniffles/ggplot/sniffles-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_ngmlr_60X.df <- read.delim(file = '4_results/subsampled/60X/sv_intersection_agreement/sniffles/ggplot/sniffles-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_pbmm2_60X.df <- read.delim(file = '4_results/subsampled/60X/sv_intersection_agreement/sniffles/ggplot/sniffles-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)

nooverlap_sniffles.df <- rbind(sniffles_minimap2_10X.df, sniffles_minimap2_20X.df, sniffles_minimap2_40X.df, sniffles_minimap2_60X.df,sniffles_ngmlr_10X.df, sniffles_ngmlr_20X.df, sniffles_ngmlr_40X.df, sniffles_ngmlr_60X.df,sniffles_pbmm2_10X.df, sniffles_pbmm2_20X.df, sniffles_pbmm2_40X.df, sniffles_pbmm2_60X.df)
col_order <- c("Depth", "Aligner", "Agreement",  "Count")
nooverlap_sniffles_reordered.df <- nooverlap_sniffles.df[, col_order]
nooverlap_sniffles_reordered.df <- nooverlap_sniffles_reordered.df[- grep("Intersection", nooverlap_sniffles_reordered.df$Agreement),]
nooverlap_sniffles_reordered.df <- nooverlap_sniffles_reordered.df[- grep("Unique", nooverlap_sniffles_reordered.df$Agreement),]
nooverlap_sniffles_reordered.df
#ggplot(nooverlap_sniffles_reordered.df, aes(x=Depth, y = Count, fill=Aligner)) + geom_col( position=position_dodge())
p_sniffles <- ggplot(nooverlap_sniffles_reordered.df, aes(x=Depth, y = Count)) + geom_col( position=position_dodge()) + facet_grid(~ Aligner, scales = "free_x", space = "free_x", switch = "x") + ylim(0, 100)
p_sniffles <- p_sniffles + xlab("Sniffles") + ylab("Non-Overlapping SV Calls (%)")
#p <- p + xlab("Sequencing Depth and Alignment Method") + ylab("Non-Overlapping SV Calls (%)") + ggtitle("Impact of Sequencing Depth on Sniffles Overlapping SV Calls") + theme(plot.title = element_text(hjust = 0.5)) #+ theme_bw()

svim_minimap2_10X.df <- read.delim(file = '4_results/subsampled/10X/sv_intersection_agreement/svim/qual_15/ggplot/svim-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_ngmlr_10X.df <- read.delim(file = '4_results/subsampled/10X/sv_intersection_agreement/svim/qual_15/ggplot/svim-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_pbmm2_10X.df <- read.delim(file = '4_results/subsampled/10X/sv_intersection_agreement/svim/qual_15/ggplot/svim-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_minimap2_20X.df <- read.delim(file = '4_results/subsampled/20X/sv_intersection_agreement/svim/qual_15/ggplot/svim-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_ngmlr_20X.df <- read.delim(file = '4_results/subsampled/20X/sv_intersection_agreement/svim/qual_15/ggplot/svim-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_pbmm2_20X.df <- read.delim(file = '4_results/subsampled/20X/sv_intersection_agreement/svim/qual_15/ggplot/svim-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_minimap2_40X.df <- read.delim(file = '4_results/subsampled/40X/sv_intersection_agreement/svim/qual_15/ggplot/svim-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_ngmlr_40X.df <- read.delim(file = '4_results/subsampled/40X/sv_intersection_agreement/svim/qual_15/ggplot/svim-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_pbmm2_40X.df <- read.delim(file = '4_results/subsampled/40X/sv_intersection_agreement/svim/qual_15/ggplot/svim-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_minimap2_60X.df <- read.delim(file = '4_results/subsampled/60X/sv_intersection_agreement/svim/qual_15/ggplot/svim-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_ngmlr_60X.df <- read.delim(file = '4_results/subsampled/60X/sv_intersection_agreement/svim/qual_15/ggplot/svim-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_pbmm2_60X.df <- read.delim(file = '4_results/subsampled/60X/sv_intersection_agreement/svim/qual_15/ggplot/svim-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)

nooverlap_svim.df <- rbind(svim_minimap2_10X.df, svim_minimap2_20X.df, svim_minimap2_40X.df, svim_minimap2_60X.df,svim_ngmlr_10X.df, svim_ngmlr_20X.df, svim_ngmlr_40X.df, svim_ngmlr_60X.df,svim_pbmm2_10X.df, svim_pbmm2_20X.df, svim_pbmm2_40X.df, svim_pbmm2_60X.df)
col_order <- c("Depth", "Aligner", "Agreement",  "Count")
nooverlap_svim_reordered.df <- nooverlap_svim.df[, col_order]
nooverlap_svim_reordered.df <- nooverlap_svim_reordered.df[- grep("Intersection", nooverlap_svim_reordered.df$Agreement),]
nooverlap_svim_reordered.df <- nooverlap_svim_reordered.df[- grep("Unique", nooverlap_svim_reordered.df$Agreement),]
#nooverlap_svim_reordered.df <- nooverlap_svim_reordered.df[- grep("Total_SVs", nooverlap_svim_reordered.df$Agreement),]

#ggplot(nooverlap_svim_reordered.df, aes(x=Depth, y = Count, fill=Aligner)) + geom_col( position=position_dodge())
p_svim <- ggplot(nooverlap_svim_reordered.df, aes(x=Depth, y = Count)) + geom_col( position=position_dodge()) + facet_grid(~ Aligner, scales = "free_x", space = "free_x", switch = "x") + ylim(0, 100)
p_svim <- p_svim + xlab("SVIM")  + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank() )
#p_svim <- p + xlab("Sequencing Depth and Alignment Method") + ylab("Non-Overlapping SV Calls (%)") + ggtitle("Impact of Sequencing Depth on SVIM Overlapping SV Calls") + theme(plot.title = element_text(hjust = 0.5)) #+ theme_bw()

pbsv_pbmm2_10X.df <- read.delim(file = '4_results/subsampled/10X/sv_intersection_agreement/pbsv/ggplot/pbsv-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
pbsv_pbmm2_20X.df <- read.delim(file = '4_results/subsampled/20X/sv_intersection_agreement/pbsv/ggplot/pbsv-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
pbsv_pbmm2_40X.df <- read.delim(file = '4_results/subsampled/40X/sv_intersection_agreement/pbsv/ggplot/pbsv-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
pbsv_pbmm2_60X.df <- read.delim(file = '4_results/subsampled/60X/sv_intersection_agreement/pbsv/ggplot/pbsv-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)

nooverlap_pbsv.df <- rbind(pbsv_pbmm2_10X.df, pbsv_pbmm2_20X.df, pbsv_pbmm2_40X.df, pbsv_pbmm2_60X.df)
col_order <- c("Depth", "Aligner", "Agreement",  "Count")
nooverlap_pbsv_reordered.df <- nooverlap_pbsv.df[, col_order]
nooverlap_pbsv_reordered.df <- nooverlap_pbsv_reordered.df[- grep("Intersection", nooverlap_pbsv_reordered.df$Agreement),]
nooverlap_pbsv_reordered.df <- nooverlap_pbsv_reordered.df[- grep("Unique", nooverlap_pbsv_reordered.df$Agreement),]
#nooverlap_pbsv_reordered.df <- nooverlap_pbsv_reordered.df[- grep("Total_SVs", nooverlap_pbsv_reordered.df$Agreement),]

p_pbsv <- ggplot(nooverlap_pbsv_reordered.df, aes(x=Depth, y = Count)) + geom_col( position=position_dodge()) + facet_grid(~ Aligner, scales = "free_x", space = "free_x", switch = "x") + ylim(0, 100)
p_pbsv <- p_pbsv + xlab("pbsv") + ylab("Non-Overlapping SV Calls (%)")

# Plot Sniffles and SVIM only
combined_plots_sniffles_svim <- p_sniffles + p_svim
combined_plots_sniffles_svim <- combined_plots_sniffles_svim + plot_annotation(title = 'Impact of Sequencing Depth on Overlapping SV Calls for Sniffles and SVIM', theme = theme(plot.title = element_text(hjust = 0.5)))
#ggarrange(p, p_svim, bp + rremove("x.text"), labels = c("A", "B"), ncol = 2, nrow = 1)

dir.create(file.path("5_plots","subsampled"), recursive = TRUE)
ggsave("3_sv_agreement_sniffles_svim_subsampled.png", plot = combined_plots_sniffles_svim, device = "png", path = "5_plots/subsampled", scale = 1, dpi = 300, limitsize = TRUE, width = 2605, height = 1715, units = "px")

# Plot pbsv, Sniffles, and SVIM

combined_plots_pbsv_sniffles_svim <- p_sniffles + p_svim + p_pbsv + plot_layout(nrow = 2, byrow = TRUE)
combined_plots_pbsv_sniffles_svim <- combined_plots_pbsv_sniffles_svim + plot_annotation(title = 'Impact of Sequencing Depth on SV Agreement', theme = theme(plot.title = element_text(hjust = 0.5)))
ggsave("3_sv_agreement_pbsv_sniffles_svim_subsampled.png", plot = combined_plots_pbsv_sniffles_svim, device = "png", path = "5_plots/subsampled", scale = 1, dpi = 300, limitsize = TRUE, width = 2605, height = 1715, units = "px")

# Plot pbsv, Sniffles, and SVIM with counts

p_sniffles_counts <- ggplot(nooverlap_sniffles_reordered.df, aes(x=Depth, y = Count, group = Depth)) + geom_col( position=position_dodge()) + facet_grid(~ Aligner, scales = "free_x", space = "free_x", switch = "x") + ylim(0, 100)
p_sniffles_counts <- p_sniffles_counts + xlab("Sniffles") + ylab("Non-Overlapping SV Calls (%)") + geom_text(aes(x=Depth, y = Count, label = Count, group = Depth), position = position_dodge(width = 1),vjust = -0.5, size = 2.5)

p_svim_counts <- ggplot(nooverlap_svim_reordered.df, aes(x=Depth, y = Count, group = Depth)) + geom_col( position=position_dodge()) + facet_grid(~ Aligner, scales = "free_x", space = "free_x", switch = "x") + ylim(0, 100)
p_svim_counts <- p_svim_counts + xlab("SVIM")  + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())  + geom_text(aes(x=Depth, y = Count, label = Count, group = Depth), position = position_dodge(width = 1),vjust = -0.5, size = 2.5)

p_pbsv_counts <- ggplot(nooverlap_pbsv_reordered.df, aes(x=Depth, y = Count, group = Depth)) + geom_col( position=position_dodge()) + facet_grid(~ Aligner, scales = "free_x", space = "free_x", switch = "x") + ylim(0, 100)
p_pbsv_counts <- p_pbsv_counts + xlab("pbsv") + ylab("Non-Overlapping SV Calls (%)") + geom_text(aes(x=Depth, y = Count, label = Count, group = Depth), position = position_dodge(width = 1),vjust = -0.5, size = 2.5)

combined_plots_pbsv_sniffles_svim_counts <- p_sniffles_counts + p_svim_counts + p_pbsv_counts + plot_layout(nrow = 2, byrow = TRUE)
combined_plots_pbsv_sniffles_svim_counts <- combined_plots_pbsv_sniffles_svim_counts + plot_annotation(title = 'Impact of Sequencing Depth on SV Agreement', theme = theme(plot.title = element_text(hjust = 0.5)))
ggsave("3_sv_agreement_pbsv_sniffles_svim_subsampled_counts.png", plot = combined_plots_pbsv_sniffles_svim_counts, device = "png", path = "5_plots/subsampled", scale = 1, dpi = 300, limitsize = TRUE, width = 2605, height = 1715, units = "px")
