# Script to recreate the results of the read1 manuscript. 

# The data for reproducing the results of the read1 manuscript are available on 
# the HPC at /rds/general/user/USERNAME/projects/neurogenomics-lab/live/Projects/TIPseq/read1_paper_data.
# The data includes paired-end and single-end BAM files for CTCF and CREB1. The 
# single-end files were derived from the paired-end files.

################################################################################
# Read 1 + 5' shift to centre cut sites
################################################################################

if(!require("remotes")) install.packages("remotes")
# remotes::install_github("neurogenomics/MotifPeeker")
library(MotifPeeker)

read1_creb_bam <- ""
read1_ctcf_bam <- ""
outdir <- ""
ctcf_motif <- universalmotif::read_jaspar("./MA1930.2.jaspar") # Motifs can be downloaded from JASPAR2024
creb_motif <- universalmotif::read_jaspar("./MA0018.5.jaspar")

MACSr::callpeak( # CREB1
  tfile = read1_creb_bam,
  nomodel = TRUE,
  qvalue = 0.01,
  shift = -75,
  extsize = 150,
  keepduplicates = "all",
  format = "BAM",
  name = "CREB_read1_shift",
  outdir = outdir
)

MACSr::callpeak( # CTCF
  tfile = read1_ctcf_bam,
  nomodel = TRUE,
  qvalue = 0.01,
  shift = -75,
  extsize = 150,
  keepduplicates = "all",
  format = "BAM",
  name = "CTCF_read1_shift",
  outdir = outdir
)

ctcf_peaks <- read_peak_file("CTCF_read1_shift_peaks.narrowPeak")
creb_peaks <- read_peak_file("CREB_read1_shift_peaks.narrowPeak")

ctcf_read1_prop <- motif_enrichment(
  peak_input = ctcf_peaks,
  motif = ctcf_motif,
  out_dir = "./",
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

creb_read1_prop <- motif_enrichment(
  peak_input = creb_peaks,
  motif = creb_motif,
  out_dir = "./",
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

ctcf_read1_sum_motif <- summit_to_motif(
  peak_input = ctcf_peaks,
  motif = ctcf_motif,
  fp_rate = 0.05,
  out_dir = "./",
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

creb_read1_sum_motif <- summit_to_motif(
  peak_input = creb_peaks,
  motif = creb_motif,
  fp_rate = 0.05,
  out_dir = "./",
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)


################################################################################
# Paired-end peak calling. format = "BAMPE". No shift. No extsize.
################################################################################

paired_creb_bam <- ""
paired_ctcf_bam <- ""

MACSr::callpeak(
  tfile = paired_creb_bam,
  qvalue = 0.01,
  keepduplicates = "all",
  format = "BAMPE", # paired-end format
  name = "CREB_paired",
  outdir = outdir
)
MACSr::callpeak(
  tfile = paired_ctcf_bam,
  qvalue = 0.01,
  keepduplicates = "all",
  format = "BAMPE",
  name = "CTCF_paired",
  outdir = outdir
)

ctcf_paired_peaks <- read_peak_file("CTCF_paired_peaks.narrowPeak")
creb_paired_peaks <- read_peak_file("CREB_paired_peaks.narrowPeak")

ctcf_paired_prop <- motif_enrichment(
  peak_input = ctcf_paired_peaks,
  motif = ctcf_motif,
  out_dir = outdir,
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

creb_paired_prop <- motif_enrichment(
  peak_input = creb_paired_peaks,
  motif = creb_motif,
  out_dir = outdir,
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

ctcf_paired_sum_motif <- summit_to_motif(
  peak_input = ctcf_paired_peaks,
  motif = ctcf_motif,
  fp_rate = 0.05,
  out_dir = outdir,
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

creb_paired_sum_motif <- summit_to_motif(
  peak_input = creb_paired_peaks,
  motif = creb_motif,
  fp_rate = 0.05,
  out_dir = outdir,
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

################################################################################
# Read 1 + NO 5' SHIFT.
################################################################################

MACSr::callpeak(
  tfile = read1_creb_bam,
  nomodel = TRUE,
  qvalue = 0.01,
  shift = 0, # leave shift at 0bp
  extsize = 150,
  keepduplicates = "all",
  format = "BAM",
  name = "CREB_noshift",
  outdir = outdir
)
MACSr::callpeak(
  tfile = read1_ctcf_bam,
  nomodel = TRUE,
  qvalue = 0.01,
  shift = 0,
  extsize = 150,
  keepduplicates = "all",
  format = "BAM",
  name = "CTCF_noshift",
  outdir = outdir
)

ctcf_read1_noshift_peaks <- read_peak_file("CTCF_noshift_peaks.narrowPeak")
creb_read1_noshift_peaks <- read_peak_file("CREB_noshift_peaks.narrowPeak")

ctcf_read1_noshift_prop <- motif_enrichment(
  peak_input = ctcf_read1_noshift_peaks,
  motif = ctcf_motif,
  out_dir = outdir,
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

creb_read1_noshift_prop <- motif_enrichment(
  peak_input = creb_read1_noshift_peaks,
  motif = creb_motif,
  out_dir = outdir,
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

ctcf_read1_noshift_sum_motif <- summit_to_motif(
  peak_input = ctcf_read1_noshift_peaks,
  motif = ctcf_motif,
  fp_rate = 0.05,
  out_dir = outdir,
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

creb_read1_noshift_sum_motif <- summit_to_motif(
  peak_input = creb_read1_noshift_peaks,
  motif = creb_motif,
  fp_rate = 0.05,
  out_dir = outdir,
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

mean(GenomicRanges::width(ctcf_peaks))
mean(GenomicRanges::width(creb_peaks))
mean(GenomicRanges::width(ctcf_read1_noshift_peaks))
mean(GenomicRanges::width(creb_read1_noshift_peaks))
mean(GenomicRanges::width(ctcf_paired_peaks))
mean(GenomicRanges::width(creb_paired_peaks))



# Plots

library(ggplot2)
library(cowplot)

ctcf_read1_sum_motif_distances <- ctcf_read1_sum_motif$distance_to_summit
creb_read1_sum_motif_distances <- creb_read1_sum_motif$distance_to_summit
ctcf_read1_noshift_sum_motif_distances <- ctcf_read1_noshift_sum_motif$distance_to_summit
creb_read1_noshift_sum_motif_distances <- creb_read1_noshift_sum_motif$distance_to_summit
ctcf_paired_sum_motif_distances <- ctcf_paired_sum_motif$distance_to_summit
creb_paired_sum_motif_distances <- creb_paired_sum_motif$distance_to_summit

# Create a data frame for each, with a condition column
ctcf_read1_sum_df <- data.frame(distance_to_summit = ctcf_read1_sum_motif_distances, condition = "CTCF first-mate [5'-shift]")
creb_read1_sum_df <- data.frame(distance_to_summit = creb_read1_sum_motif_distances, condition = "CREB first-mate [5'-shift]")
ctcf_read1_noshift_df <- data.frame(distance_to_summit = ctcf_read1_noshift_sum_motif_distances, condition = "CTCF first-mate [No 5'-shift]")
creb_read1_noshift_df <- data.frame(distance_to_summit = creb_read1_noshift_sum_motif_distances, condition = "CREB first-mate [No 5'-shift]")
ctcf_paired_df <- data.frame(distance_to_summit = ctcf_paired_sum_motif_distances, condition = "CTCF paired-end")
creb_paired_df <- data.frame(distance_to_summit = creb_paired_sum_motif_distances, condition = "CREB paired-end")

combined_df <- rbind(ctcf_read1_sum_df, creb_read1_sum_df, ctcf_read1_noshift_df, creb_read1_noshift_df, ctcf_paired_df, creb_paired_df)
combined_df$condition <- factor(combined_df$condition, levels = c(
  "CTCF first-mate [5'-shift]", "CREB first-mate [5'-shift]",
  "CTCF first-mate [No 5'-shift]", "CREB first-mate [No 5'-shift]",
  "CTCF paired-end", "CREB paired-end"
))
ggplot(combined_df, aes(x = distance_to_summit)) +
  geom_line(stat = "density", linetype = "solid", linewidth = 0.8) +
  facet_wrap(~ condition, scales = "free_y", ncol = 2) + # Facets arranged by 'condition'
  labs(
    x = "Summit-to-motif distance (bp)",
    y = "Density",
  ) +
  theme_minimal() +
  theme(
    plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
    axis.title.x = element_text(size = 16, margin = margin(t = 10), colour = "black"),
    axis.title.y = element_text(size = 16, margin = margin(r = 10), colour = "black"),
    axis.text = element_text(size = 13, colour = "black"),
    strip.text = element_text(size = 14, face = "bold", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.5)
  ) +
  expand_limits(x = 0, y = -1e-04) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(-400,400),
                     breaks = seq(-1e+05, 1e+05, by = 100)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

# The de novo motif discovery results are on the HPC at 
# /rds/general/user/USERNAME/projects/neurogenomics-lab/live/Projects/TIPseq/read1_paper_data
# library(MotifStats)
# 
# ctcf_peaks <- read_peak_file("./CTCF_peaks.narrowPeak")
# genome_build <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
# ctcf_sequences <- BSgenome::getSeq(genome_build,
#                                    ctcf_peaks)
# 
# # Function to adjust ranges within sequence boundaries
# adjust_ranges_within_bounds <- function(peaks, genome_build) {
#   adjusted_ranges <- IRanges(
#     start = pmax(peaks$summit - 25, 1),  # Ensure start is not less than 1. 
#     end = pmin(peaks$summit + 25, GenomeInfoDb::seqlengths(genome_build)[as.character(seqnames(ctcf_peaks))])  # Ensure end is not more than the sequence length
#   )
#   return(adjusted_ranges)
# }
# 
# new_ranges <- adjust_ranges_within_bounds(ctcf_peaks, genome_build)
# ranges(ctcf_peaks) <- new_ranges
# ctcf_sequences <- BSgenome::getSeq(genome_build,
#                                    ctcf_peaks)
# 
# run_time <- Sys.time()
# memes::runStreme(input = ctcf_sequences,
#                  control = "shuffle", # preserve trinucleotide frequencies in the shuffled sequences.
#                  outdir = "./ctcf_read1_shift_25plusminus",
#                  nmotifs = 5, # max number of motifs to discover
#                  minw = 16, # min width of a motif
#                  maxw = 24) # max width of a motif
# # The CTCF motif is quite long.
# run_time <- Sys.time() - run_time
