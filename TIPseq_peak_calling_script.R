################################################################################
# Read 1 + 5' shift to centre cut sites
################################################################################

# if(!require("remotes")) install.packages("remotes")
# remotes::install_github("neurogenomics/MotifStats")
library(MotifStats)

MACSr::callpeak( # CREB1
  tfile = "./read1_creb.bam",
  nomodel = TRUE,
  qvalue = 0.01,
  shift = -75,
  extsize = 150,
  keepduplicates = "all",
  format = "BAM",
  name = "CREB_read1_shift",
  outdir = "./"
)

MACSr::callpeak( # CTCF
  tfile = "./read1_ctcf.bam",
  nomodel = TRUE,
  qvalue = 0.01,
  shift = -75,
  extsize = 150,
  keepduplicates = "all",
  format = "BAM",
  name = "CTCF_read1_shift",
  outdir = "./"
)

ctcf_peaks <- read_peak_file("./CTCF_read1_shift_peaks.narrowPeak")
creb_peaks <- read_peak_file("./CREB_read1_shift_peaks.narrowPeak")
ctcf_motif <- universalmotif::read_jaspar("./MA1930.2.jaspar") # Motifs can be downloaded from JASPAR2024
creb_motif <- universalmotif::read_jaspar("./MA0018.5.jaspar")

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
# Paired-end peak calling. format = "BAMPE".
################################################################################

MACSr::callpeak(
  tfile = "./SK032_1_SK035_1_SK035_11_R1.target.linear_dedup.sorted.bam",
  qvalue = 0.01,
  keepduplicates = "all",
  format = "BAMPE", # paired-end format
  name = "CREB_paired",
  outdir = "./"
)
MACSr::callpeak(
  tfile = "./CTCF_R1.target.linear_dedup.sorted.bam",
  qvalue = 0.01,
  keepduplicates = "all",
  format = "BAMPE",
  name = "CTCF_paired",
  outdir = "./"
)

ctcf_paired_peaks <- read_peak_file("./CTCF_paired_peaks.narrowPeak")
creb_paired_peaks <- read_peak_file("./CREB_paired_peaks.narrowPeak")

ctcf_paired_prop <- motif_enrichment(
  peak_input = ctcf_paired_peaks,
  motif = ctcf_motif,
  out_dir = ".",
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

creb_paired_prop <- motif_enrichment(
  peak_input = creb_paired_peaks,
  motif = creb_motif,
  out_dir = ".",
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

ctcf_paired_sum_motif <- summit_to_motif(
  peak_input = ctcf_paired_peaks,
  motif = ctcf_motif,
  fp_rate = 0.05,
  out_dir = ".",
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

creb_paired_sum_motif <- summit_to_motif(
  peak_input = creb_paired_peaks,
  motif = creb_motif,
  fp_rate = 0.05,
  out_dir = ".",
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

################################################################################
# Read 1 + NO 5' SHIFT.
################################################################################

MACSr::callpeak(
  tfile = "./read1_creb.bam",
  nomodel = TRUE,
  qvalue = 0.01,
  shift = 0, # leave shift at 0bp
  extsize = 150,
  keepduplicates = "all",
  format = "BAM",
  name = "CREB_noshift",
  outdir = "./"
)
MACSr::callpeak(
  tfile = "./read1_ctcf.bam",
  nomodel = TRUE,
  qvalue = 0.01,
  shift = 0,
  extsize = 150,
  keepduplicates = "all",
  format = "BAM",
  name = "CTCF_noshift",
  outdir = "./"
)

ctcf_read1_noshift_peaks <- read_peak_file("./CTCF_noshift_peaks.narrowPeak")
creb_read1_noshift_peaks <- read_peak_file("./CREB_noshift_peaks.narrowPeak")

ctcf_read1_noshift_prop <- motif_enrichment(
  peak_input = ctcf_read1_noshift_peaks,
  motif = ctcf_motif,
  out_dir = ".",
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

creb_read1_noshift_prop <- motif_enrichment(
  peak_input = creb_read1_noshift_peaks,
  motif = creb_motif,
  out_dir = ".",
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

ctcf_read1_noshift_sum_motif <- summit_to_motif(
  peak_input = ctcf_read1_noshift_peaks,
  motif = ctcf_motif,
  fp_rate = 0.05,
  out_dir = ".",
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

creb_read1_noshift_sum_motif <- summit_to_motif(
  peak_input = creb_read1_noshift_peaks,
  motif = creb_motif,
  fp_rate = 0.05,
  out_dir = ".",
  genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

mean(GenomicRanges::width(ctcf_peaks))
mean(GenomicRanges::width(creb_peaks))
mean(GenomicRanges::width(ctcf_read1_noshift_peaks))
mean(GenomicRanges::width(creb_read1_noshift_peaks))
mean(GenomicRanges::width(ctcf_paired_peaks))
mean(GenomicRanges::width(creb_paired_peaks))