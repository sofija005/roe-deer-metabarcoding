library(dada2)

base <- "/home/sofija/INTERNSHIP/samples/Capreolus_capreolus/Fungi_ITS"

filt_path <- file.path(base, "qualFiltered_out")
error_path <- file.path(base, "dada2_out", "error_models")
intermediate_path <- file.path(base, "dada2_out", "intermediate")
raw_path <- file.path(base, "dada2_out", "raw")
tracking_path <- file.path(base, "dada2_out", "tracking")

dir.create(error_path, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_path, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_path, recursive = TRUE, showWarnings = FALSE)
dir.create(tracking_path, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Locate quality-filtered reads
# ------------------------------------------------------------

fnFs <- sort(list.files(
  filt_path,
  pattern = "_R1\\.fq\\.gz$",
  full.names = TRUE
))

fnRs <- sort(list.files(
  filt_path,
  pattern = "_R2\\.fq\\.gz$",
  full.names = TRUE
))

stopifnot(length(fnFs) == length(fnRs))
stopifnot(length(fnFs) == 20)

sample.names <- sub("_R1\\.fq\\.gz$", "", basename(fnFs))

names(fnFs) <- sample.names
names(fnRs) <- sample.names

cat("Samples detected:", length(sample.names), "\n")
print(sample.names)

# ------------------------------------------------------------
# Learn error rates
# ------------------------------------------------------------

cat("\nLearning forward error rates...\n")

errF <- learnErrors(
  fnFs,
  multithread = TRUE,
  randomize = TRUE
)

cat("\nLearning reverse error rates...\n")

errR <- learnErrors(
  fnRs,
  multithread = TRUE,
  randomize = TRUE
)

saveRDS(
  errF,
  file.path(intermediate_path, "LIB3_errF.rds")
)

saveRDS(
  errR,
  file.path(intermediate_path, "LIB3_errR.rds")
)

# Error plots

pdf(
  file.path(error_path, "LIB3_error_rates_R1.pdf"),
  width = 12,
  height = 9
)
print(plotErrors(errF, nominalQ = TRUE))
dev.off()

pdf(
  file.path(error_path, "LIB3_error_rates_R2.pdf"),
  width = 12,
  height = 9
)
print(plotErrors(errR, nominalQ = TRUE))
dev.off()

# ------------------------------------------------------------
# Denoising
# ------------------------------------------------------------

cat("\nDenoising R1...\n")

dadaFs <- dada(
  fnFs,
  err = errF,
  multithread = TRUE
)

cat("\nDenoising R2...\n")

dadaRs <- dada(
  fnRs,
  err = errR,
  multithread = TRUE
)

saveRDS(
  dadaFs,
  file.path(intermediate_path, "LIB3_dadaFs.rds")
)

saveRDS(
  dadaRs,
  file.path(intermediate_path, "LIB3_dadaRs.rds")
)

# ------------------------------------------------------------
# Merge paired-end reads
# ------------------------------------------------------------

cat("\nMerging paired-end reads...\n")

mergers <- mergePairs(
  dadaFs,
  fnFs,
  dadaRs,
  fnRs,
  verbose = TRUE
)

saveRDS(
  mergers,
  file.path(intermediate_path, "LIB3_mergers.rds")
)

# ------------------------------------------------------------
# Sequence table
# ------------------------------------------------------------

seqtab <- makeSequenceTable(mergers)

saveRDS(
  seqtab,
  file.path(raw_path, "LIB3_seqtab_raw.rds")
)

# ------------------------------------------------------------
# Assign stable ASV IDs
# ------------------------------------------------------------

asv_sequences <- colnames(seqtab)

asv_ids <- paste0(
  "ASV",
  seq_along(asv_sequences)
)

raw_map <- data.frame(
  ASV_ID = asv_ids,
  Sequence = asv_sequences
)

write.csv(
  raw_map,
  file.path(raw_path, "LIB3_ASV_sequences_raw.csv"),
  row.names = FALSE
)

# ASV abundance table: ASVs x samples

asv_table <- as.data.frame(t(seqtab))
rownames(asv_table) <- asv_ids

write.csv(
  asv_table,
  file.path(raw_path, "LIB3_ASV_table_raw.csv"),
  row.names = TRUE
)

# Raw ASV FASTA

fasta_lines <- as.vector(
  rbind(
    paste0(">", asv_ids),
    asv_sequences
  )
)

writeLines(
  fasta_lines,
  file.path(raw_path, "LIB3_ASVs_raw.fasta")
)

# ------------------------------------------------------------
# Tracking
# ------------------------------------------------------------

getN <- function(x) sum(getUniques(x))

denoised_F <- sapply(dadaFs, getN)
denoised_R <- sapply(dadaRs, getN)
merged <- sapply(mergers, function(x) sum(x$abundance))

filter_summary <- read.csv(
  file.path(filt_path, "seq_count_summary.csv"),
  stringsAsFactors = FALSE
)

tracking <- data.frame(
  Sample = sample.names,
  Input = filter_summary$reads.in[
    match(sample.names, filter_summary$sample)
  ],
  Quality_filtered = filter_summary$reads.out[
    match(sample.names, filter_summary$sample)
  ],
  Denoised_R1 = denoised_F[sample.names],
  Denoised_R2 = denoised_R[sample.names],
  Merged = merged[sample.names]
)

tracking$Merge_retention_percent <- round(
  100 * tracking$Merged /
    pmin(tracking$Denoised_R1, tracking$Denoised_R2),
  2
)

write.csv(
  tracking,
  file.path(tracking_path, "LIB3_dada2_part1_tracking.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("LIB3 DADA2 PART 1 COMPLETED\n")
cat("============================================\n")

print(tracking)

cat("\nTOTALS\n")
cat("Input:", sum(tracking$Input), "\n")
cat("Quality-filtered:", sum(tracking$Quality_filtered), "\n")
cat("Denoised R1:", sum(tracking$Denoised_R1), "\n")
cat("Denoised R2:", sum(tracking$Denoised_R2), "\n")
cat("Merged:", sum(tracking$Merged), "\n")

cat(
  "Overall merge retention:",
  round(
    100 * sum(tracking$Merged) /
      min(
        sum(tracking$Denoised_R1),
        sum(tracking$Denoised_R2)
      ),
    2
  ),
  "%\n"
)

cat("Raw ASVs:", ncol(seqtab), "\n")

cat("\nSTOP AND REVIEW BEFORE CHIMERA REMOVAL.\n")

