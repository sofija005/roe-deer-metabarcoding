library(dada2)

base <- "/home/sofija/INTERNSHIP/samples/Capreolus_capreolus/Plants_ITS"

# Load objects
dadaFs <- readRDS(file.path(
  base, "dada2_out/intermediate/LIB1_dadaFs.rds"
))

dadaRs <- readRDS(file.path(
  base, "dada2_out/intermediate/LIB1_dadaRs.rds"
))

mergers <- readRDS(file.path(
  base, "dada2_out/intermediate/LIB1_mergers.rds"
))

seqtab <- readRDS(file.path(
  base, "dada2_out/raw/LIB1_seqtab_raw.rds"
))

seqtab.nochim <- readRDS(file.path(
  base, "dada2_out/nochim/LIB1_seqtab_nochim.rds"
))

seqtab.uncross <- readRDS(file.path(
  base, "dada2_out/uncross2/LIB1_seqtab_uncross2.rds"
))

# Filtering summary
f <- read.csv(
  file.path(base, "qualFiltered_out/seq_count_summary.csv"),
  check.names = FALSE
)

names(f)[1] <- "file"
f$Sample <- sub("_R1\\.fq\\.gz$", "", f$file)

samples <- names(dadaFs)

getN <- function(x) sum(getUniques(x))

# Per-sample tracking
tracking <- data.frame(
  Sample = samples,
  Input = f$reads.in[match(samples, f$Sample)],
  Quality_filtered = f$reads.out[match(samples, f$Sample)],
  Denoised_R1 = sapply(dadaFs, getN),
  Denoised_R2 = sapply(dadaRs, getN),
  Merged = rowSums(seqtab),
  Non_chimeric = rowSums(seqtab.nochim),
  UNCROSS2_filtered = rowSums(seqtab.uncross)
)

tracking$Contaminant_filtered <- NA

write.csv(
  tracking,
  file.path(base, "dada2_out/tracking/LIB1_final_sequence_tracking.csv"),
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
  file.path(base, "dada2_out/tracking/LIB1_ASV_tracking.csv"),
  row.names = FALSE
)

cat("\n=== FINAL LIB1 SEQUENCE TRACKING ===\n")
print(tracking)

cat("\n=== FINAL LIB1 ASV TRACKING ===\n")
print(asv_tracking)

cat("\nContaminant assessment: NOT PERFORMED - no suitable controls/metadata.\n")
