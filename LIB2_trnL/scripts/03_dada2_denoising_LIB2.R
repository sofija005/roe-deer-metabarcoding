library(dada2)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

base <- "/home/sofija/INTERNSHIP/samples/Capreolus_capreolus/Plants_trnLc-trnLh"

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
# Find quality-filtered paired-end reads
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

sampleF <- sub("_R1\\.fq\\.gz$", "", basename(fnFs))
sampleR <- sub("_R2\\.fq\\.gz$", "", basename(fnRs))

stopifnot(length(fnFs) == length(fnRs))
stopifnot(identical(sampleF, sampleR))

samples <- sampleF

cat("Number of samples:", length(samples), "\n")
cat("Samples:\n")
print(samples)


# ------------------------------------------------------------
# 1. Learn error rates
# ------------------------------------------------------------

cat("\nLearning forward-read error rates...\n")

errF <- learnErrors(
  fnFs,
  multithread = TRUE
)

cat("\nLearning reverse-read error rates...\n")

errR <- learnErrors(
  fnRs,
  multithread = TRUE
)

saveRDS(
  errF,
  file.path(error_path, "LIB2_errF.rds")
)

saveRDS(
  errR,
  file.path(error_path, "LIB2_errR.rds")
)


# ------------------------------------------------------------
# Error-rate plots
# ------------------------------------------------------------

pdf(file.path(error_path, "LIB2_error_rates_R1.pdf"))
plotErrors(errF, nominalQ = TRUE)
dev.off()

pdf(file.path(error_path, "LIB2_error_rates_R2.pdf"))
plotErrors(errR, nominalQ = TRUE)
dev.off()


# ------------------------------------------------------------
# 2. Denoising
# ------------------------------------------------------------

cat("\nDenoising forward reads...\n")

dadaFs <- dada(
  fnFs,
  err = errF,
  multithread = TRUE
)

cat("\nDenoising reverse reads...\n")

dadaRs <- dada(
  fnRs,
  err = errR,
  multithread = TRUE
)

names(dadaFs) <- samples
names(dadaRs) <- samples

saveRDS(
  dadaFs,
  file.path(intermediate_path, "LIB2_dadaFs.rds")
)

saveRDS(
  dadaRs,
  file.path(intermediate_path, "LIB2_dadaRs.rds")
)


# ------------------------------------------------------------
# 3. Merge paired-end reads
# ------------------------------------------------------------

cat("\nMerging paired-end reads...\n")

mergers <- mergePairs(
  dadaFs,
  fnFs,
  dadaRs,
  fnRs,
  verbose = TRUE
)

names(mergers) <- samples

saveRDS(
  mergers,
  file.path(intermediate_path, "LIB2_mergers.rds")
)


# ------------------------------------------------------------
# 4. Construct RAW ASV sequence table
# ------------------------------------------------------------

seqtab <- makeSequenceTable(mergers)

saveRDS(
  seqtab,
  file.path(raw_path, "LIB2_seqtab_raw.rds")
)


# ------------------------------------------------------------
# Stable ASV identifiers
# ------------------------------------------------------------

asv_sequences <- colnames(seqtab)
asv_ids <- paste0("ASV", seq_along(asv_sequences))

asv_map <- data.frame(
  ASV_ID = asv_ids,
  Sequence = asv_sequences
)

write.csv(
  asv_map,
  file.path(raw_path, "LIB2_ASV_sequences_raw.csv"),
  row.names = FALSE
)


# ASV abundance table
asv_table <- as.data.frame(t(seqtab))
rownames(asv_table) <- asv_ids

write.csv(
  asv_table,
  file.path(raw_path, "LIB2_ASV_table_raw.csv"),
  row.names = TRUE
)


# RAW ASV FASTA
fasta_lines <- as.vector(
  rbind(
    paste0(">", asv_ids),
    asv_sequences
  )
)

writeLines(
  fasta_lines,
  file.path(raw_path, "LIB2_ASVs_raw.fasta")
)


# ------------------------------------------------------------
# Tracking
# ------------------------------------------------------------

getN <- function(x) sum(getUniques(x))

filter_summary <- read.csv(
  file.path(filt_path, "seq_count_summary.csv"),
  stringsAsFactors = FALSE
)

stopifnot(all(samples %in% filter_summary$sample))

tracking <- data.frame(
  Sample = samples,
  Input = filter_summary$reads.in[
    match(samples, filter_summary$sample)
  ],
  Quality_filtered = filter_summary$reads.out[
    match(samples, filter_summary$sample)
  ],
  Denoised_R1 = sapply(dadaFs, getN),
  Denoised_R2 = sapply(dadaRs, getN),
  Merged = sapply(mergers, getN)
)

tracking$Merge_retained_percent <- round(
  100 * tracking$Merged /
    pmin(tracking$Denoised_R1, tracking$Denoised_R2),
  2
)

write.csv(
  tracking,
  file.path(tracking_path, "LIB2_tracking_part1.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("LIB2 DADA2 PART 1 COMPLETED\n")
cat("============================================\n")

cat("\nRaw ASVs:", ncol(seqtab), "\n")

cat("\nRead tracking:\n")
print(tracking)

cat("\nTOTALS:\n")
cat("Input:", sum(tracking$Input), "\n")
cat("Quality-filtered:", sum(tracking$Quality_filtered), "\n")
cat("Denoised R1:", sum(tracking$Denoised_R1), "\n")
cat("Denoised R2:", sum(tracking$Denoised_R2), "\n")
cat("Merged:", sum(tracking$Merged), "\n")

cat(
  "Overall merging retention:",
  round(
    100 * sum(tracking$Merged) /
      min(sum(tracking$Denoised_R1),
          sum(tracking$Denoised_R2)),
    2
  ),
  "%\n"
)

cat("\nSTOP HERE AND REVIEW MERGING RESULTS BEFORE CHIMERA REMOVAL.\n")

