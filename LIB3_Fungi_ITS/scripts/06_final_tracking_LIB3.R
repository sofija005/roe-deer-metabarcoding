library(dada2)

base <- "/home/sofija/INTERNSHIP/samples/Capreolus_capreolus/Fungi_ITS"

tracking_path <- file.path(base, "dada2_out", "tracking")
dir.create(tracking_path, recursive = TRUE, showWarnings = FALSE)

# Load DADA2 outputs
dadaFs <- readRDS(file.path(
  base, "dada2_out/intermediate/LIB3_dadaFs.rds"
))

dadaRs <- readRDS(file.path(
  base, "dada2_out/intermediate/LIB3_dadaRs.rds"
))

seqtab <- readRDS(file.path(
  base, "dada2_out/raw/LIB3_seqtab_raw.rds"
))

seqtab.nochim <- readRDS(file.path(
  base, "dada2_out/nochim/LIB3_seqtab_nochim.rds"
))

seqtab.uncross <- readRDS(file.path(
  base, "dada2_out/uncross2/LIB3_seqtab_uncross2.rds"
))

# Filtering summary
f <- read.csv(
  file.path(base, "qualFiltered_out/seq_count_summary.csv"),
  stringsAsFactors = FALSE
)

samples <- names(dadaFs)

getN <- function(x) sum(getUniques(x))

# Final sequence tracking
tracking <- data.frame(
  Sample = samples,

  Input = f$reads.in[
    match(samples, f$sample)
  ],

  Quality_filtered = f$reads.out[
    match(samples, f$sample)
  ],

  Denoised_R1 = sapply(dadaFs, getN),

  Denoised_R2 = sapply(dadaRs, getN),

  Merged = rowSums(seqtab)[samples],

  Non_chimeric = rowSums(seqtab.nochim)[samples],

  UNCROSS2_filtered = rowSums(seqtab.uncross)[samples]
)

# No contaminant filtering
tracking$Contaminant_filtered <- NA

write.csv(
  tracking,
  file.path(tracking_path, "LIB3_final_sequence_tracking.csv"),
  row.names = FALSE
)

# ASV tracking
asv_tracking <- data.frame(
  Step = c(
    "Raw ASV table",
    "Non-chimeric",
    "UNCROSS2-filtered",
    "Final"
  ),

  ASVs = c(
    ncol(seqtab),
    ncol(seqtab.nochim),
    ncol(seqtab.uncross),
    ncol(seqtab.uncross)
  )
)

write.csv(
  asv_tracking,
  file.path(tracking_path, "LIB3_ASV_tracking.csv"),
  row.names = FALSE
)

# Print results
cat("\n=== FINAL LIB3 SEQUENCE TRACKING ===\n")
print(tracking)

cat("\n=== FINAL LIB3 ASV TRACKING ===\n")
print(asv_tracking)

cat("\n=== TOTAL READS ===\n")
cat("Input:", sum(tracking$Input), "\n")
cat("Quality-filtered:", sum(tracking$Quality_filtered), "\n")
cat("Denoised R1:", sum(tracking$Denoised_R1), "\n")
cat("Denoised R2:", sum(tracking$Denoised_R2), "\n")
cat("Merged:", sum(tracking$Merged), "\n")
cat("Non-chimeric:", sum(tracking$Non_chimeric), "\n")
cat("UNCROSS2-filtered:", sum(tracking$UNCROSS2_filtered), "\n")

cat(
  "\nContaminant assessment: NOT PERFORMED - no suitable controls/metadata.\n"
)

