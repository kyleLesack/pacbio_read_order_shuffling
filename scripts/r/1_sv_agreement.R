library(ggplot2)
library(tidyverse)
setwd("/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling")

pbsv_pbmm2 <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/pbsv/ggplot/pbsv-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_minimap2 <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/sniffles/ggplot/sniffles-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_ngmlr <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/sniffles/ggplot/sniffles-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_pbmm2 <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/sniffles/ggplot/sniffles-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_minimap2 <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/svim/qual_15/ggplot/svim-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_ngmlr <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/svim/qual_15/ggplot/svim-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_pbmm2 <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/svim/qual_15/ggplot/svim-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)

pbsv_pbmm2.df <- data.frame(pbsv_pbmm2)
sniffles_minimap2.df <- data.frame(sniffles_minimap2)
sniffles_ngmlr.df <- data.frame(sniffles_ngmlr)
sniffles_pbmm2.df <- data.frame(sniffles_pbmm2)
svim_minimap2.df <- data.frame(svim_minimap2)
svim_ngmlr.df <- data.frame(svim_ngmlr)
svim_pbmm2.df <- data.frame(svim_pbmm2)

intersection.df <- rbind(pbsv_pbmm2.df, sniffles_minimap2.df, sniffles_ngmlr.df, sniffles_pbmm2.df,svim_minimap2.df, svim_ngmlr.df, svim_pbmm2.df)
col_order <- c("Aligner", "Caller", "Agreement", "Count")
intersection_reordered.df <- intersection.df[, col_order]
intersection_reordered.df <- intersection_reordered.df[- grep("Non-Overlapping", intersection_reordered.df$Agreement),]
intersection_reordered.df <- intersection_reordered.df[- grep("Total_SVs", intersection_reordered.df$Agreement),]
intersection_reordered.df <- intersection_reordered.df %>% mutate_all(funs(str_replace(., "Intersection", "Intersecting calls")))
intersection_reordered.df <- intersection_reordered.df %>% mutate_all(funs(str_replace(., "Unique", "Unique calls")))
intersection_reordered.df <- intersection_reordered.df %>% mutate_all(funs(str_replace(., "sniffles", "Sniffles")))
intersection_reordered.df <- intersection_reordered.df %>% mutate_all(funs(str_replace(., "svim", "SVIM")))
intersection_reordered.df <- intersection_reordered.df %>% mutate_all(funs(str_replace(., "minimap2", "Minimap2")))
intersection_reordered.df <- intersection_reordered.df %>% mutate_all(funs(str_replace(., "ngmlr", "NGMLR")))
intersection_reordered.df$Count = as.numeric(intersection_reordered.df$Count)

p <- ggplot(intersection_reordered.df, aes(x=Caller, y = Count, fill=Agreement)) + geom_col( position=position_dodge()) + facet_grid(~ Aligner, scales = "free_x", space = "free_x", switch = "x")
p <- p + xlab("SV Calling and Alignment Method") + ggtitle("Intersection of Predicted Structural Variants") + theme(plot.title = element_text(hjust = 0.5)) #+ theme_bw()

dir.create(file.path("5_plots","full_depth"), recursive = TRUE)
ggsave("1_sv_agreement.png", plot = p, device = "png", path = "5_plots/full_depth", scale = 1, dpi = 300, limitsize = TRUE, width = 2605, height = 1715, units = "px")

# Plot with values included on top of bars
p2 <- ggplot(intersection_reordered.df, aes(x=Caller, y = Count, fill=Agreement, group = Agreement)) + geom_col( position=position_dodge()) + facet_grid(~ Aligner, scales = "free_x", space = "free_x", switch = "x")
p2 <- p2 + xlab("SV Calling and Alignment Method") + ggtitle("Intersection of Predicted Structural Variants") + theme(plot.title = element_text(hjust = 0.5)) + geom_text(aes(x=Caller, y = Count, label = Count, group = Agreement), position = position_dodge(width = 1),vjust = -0.5, size = 2.5)
ggsave("1_sv_agreement-counts.png", plot = p2, device = "png", path = "5_plots/full_depth", scale = 1, dpi = 300, limitsize = TRUE, width = 2605, height = 1715, units = "px")


