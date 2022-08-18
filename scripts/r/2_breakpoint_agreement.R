library(ggplot2)
library(tidyverse)
setwd("/bulk/worm_lab/mrkyle/pacbio_read_order_shuffling")

pbsv_pbmm2.df <- read.delim(file = '4_results/full_depth/breakpoint_agreement/pbsv/ggplot/pbsv-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_minimap2.df <- read.delim(file = '4_results/full_depth/breakpoint_agreement/sniffles/ggplot/sniffles-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_ngmlr.df <- read.delim(file = '4_results/full_depth/breakpoint_agreement/sniffles/ggplot/sniffles-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_pbmm2.df <- read.delim(file = '4_results/full_depth/breakpoint_agreement/sniffles/ggplot/sniffles-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_minimap2.df <- read.delim(file = '4_results/full_depth/breakpoint_agreement/svim/qual_15/ggplot/svim-minimap2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_ngmlr.df <- read.delim(file = '4_results/full_depth/breakpoint_agreement/svim/qual_15/ggplot/svim-ngmlr.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_pbmm2.df <- read.delim(file = '4_results/full_depth/breakpoint_agreement/svim/qual_15/ggplot/svim-pbmm2.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)

breakpoint_agreement.df <- rbind(pbsv_pbmm2.df, sniffles_minimap2.df, sniffles_ngmlr.df, sniffles_pbmm2.df,svim_minimap2.df, svim_ngmlr.df, svim_pbmm2.df)
col_order <- c("Aligner", "Caller", "Agreement", "Count")
breakpoint_agreement_reordered.df <- breakpoint_agreement.df[, col_order]
breakpoint_agreement_reordered.df <- breakpoint_agreement_reordered.df %>% mutate_all(funs(str_replace(., "sniffles", "Sniffles")))
breakpoint_agreement_reordered.df <- breakpoint_agreement_reordered.df %>% mutate_all(funs(str_replace(., "svim", "SVIM")))
breakpoint_agreement_reordered.df <- breakpoint_agreement_reordered.df %>% mutate_all(funs(str_replace(., "minimap2", "Minimap2")))
breakpoint_agreement_reordered.df <- breakpoint_agreement_reordered.df %>% mutate_all(funs(str_replace(., "ngmlr", "NGMLR")))
breakpoint_agreement_reordered.df <- breakpoint_agreement_reordered.df[- grep("Discordant", breakpoint_agreement_reordered.df$Agreement),]
breakpoint_agreement_reordered.df$Count = as.numeric(breakpoint_agreement_reordered.df$Count)

# Convert Agreement column to a factor and reorder so that Same breakpoints occurs first.
# Reordering makes Same breakpoints bar to occur first in plot. 
# This is more consistent with the SV agreement plot, which shows the intersection count before the unique count.
breakpoint_agreement_reordered.df$Agreement <- as.factor(breakpoint_agreement_reordered.df$Agreement)
breakpoint_agreement_reordered.df$Agreement <- factor(breakpoint_agreement_reordered.df$Agreement, levels = c("Same breakpoints", "Different breakpoints"))

p <- ggplot(breakpoint_agreement_reordered.df, aes(x=Caller, y = Count, fill=Agreement)) + geom_col( position=position_dodge()) + facet_grid(~ Aligner, scales = "free_x", space = "free_x", switch = "x")
p <- p + xlab("SV Calling and Alignment Method") + ggtitle("Breakpoint Concordance for Overlapping Predictions ") + theme(plot.title = element_text(hjust = 0.5)) #+ theme_bw()

dir.create(file.path("5_plots","full_depth"), recursive = TRUE)
ggsave("2_breakpoint_agreement.png", plot = p, device = "png", path = "5_plots/full_depth", scale = 1, dpi = 300, limitsize = TRUE, width = 2605, height = 1715, units = "px")

# Plots with values included on top of bars

p2 <- ggplot(breakpoint_agreement_reordered.df, aes(x=Caller, y = Count, fill=Agreement, group = Agreement)) + geom_col( position=position_dodge()) + facet_grid(~ Aligner, scales = "free_x", space = "free_x", switch = "x")
p2 <- p2 + xlab("SV Calling and Alignment Method") + ggtitle("Breakpoint Concordance for Overlapping Predictions ") + theme(plot.title = element_text(hjust = 0.5)) + geom_text(aes(x=Caller, y = Count, label = Count, group = Agreement), position = position_dodge(width = 1),vjust = -0.5, size = 2.5)
ggsave("2_breakpoint_agreement-counts.png", plot = p2, device = "png", path = "5_plots/full_depth", scale = 1, dpi = 300, limitsize = TRUE, width = 2605, height = 1715, units = "px")

