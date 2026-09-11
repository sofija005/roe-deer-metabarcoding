library(data.table)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

base <- "/home/sofija/INTERNSHIP/samples/Capreolus_capreolus/Plants_trnLc-trnLh"

raw_path <- file.path(base, "dada2_out", "raw")
nochim_path <- file.path(base, "dada2_out", "nochim")
out_path <- file.path(base, "dada2_out", "uncross2")

dir.create(out_path, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# Parameters from Bioscanflow
# ------------------------------------------------------------

set_f <- 0.03
set_p <- 1


# ------------------------------------------------------------
# Load non-chimeric ASV table
# ------------------------------------------------------------

seqtab.nochim <- readRDS(
  file.path(nochim_path, "LIB2_seqtab_nochim.rds")
)

raw_map <- fread(
  file.path(raw_path, "LIB2_ASV_sequences_raw.csv")
)


# ------------------------------------------------------------
# Match sequences to original ASV IDs
# ------------------------------------------------------------

idx <- match(colnames(seqtab.nochim), raw_map$Sequence)

stopifnot(!any(is.na(idx)))

asv_ids <- raw_map$ASV_ID[idx]
asv_sequences <- colnames(seqtab.nochim)


# ------------------------------------------------------------
# Convert to ASV x sample table
# ------------------------------------------------------------

wide <- as.data.table(t(seqtab.nochim))

wide[, ASV := asv_ids]

setcolorder(
  wide,
  c("ASV", setdiff(names(wide), "ASV"))
)


# ------------------------------------------------------------
# Convert to long format
# ------------------------------------------------------------

ASVTAB <- melt(
  wide,
  id.vars = "ASV",
  variable.name = "SampleID",
  value.name = "Abundance"
)

ASVTAB <- ASVTAB[Abundance > 0]


# ------------------------------------------------------------
# Total abundance per ASV
# ------------------------------------------------------------

ASVTAB[, Total := sum(Abundance, na.rm = TRUE), by = ASV]


# ------------------------------------------------------------
# UNCROSS2 scoring function
# ------------------------------------------------------------

uncross_score <- function(x, N, n, f = 0.01, tmin = 0.1, p = 1) {

  z <- f * N / n

  sc <- 2 / (1 + exp(x / z)^p)

  data.table(
    Score = sc,
    TagJump = sc >= tmin
  )
}


# ------------------------------------------------------------
# Calculate UNCROSS2 scores
# ------------------------------------------------------------

n_samples <- uniqueN(ASVTAB$SampleID)

scores <- uncross_score(
  x = ASVTAB$Abundance,
  N = ASVTAB$Total,
  n = n_samples,
  f = set_f,
  p = set_p
)

ASVTAB <- cbind(ASVTAB, scores)


# Save complete score table
fwrite(
  ASVTAB,
  file.path(out_path, "LIB2_UNCROSS2_scores.tsv"),
  sep = "\t"
)


# ------------------------------------------------------------
# Statistics before filtering
# ------------------------------------------------------------

asvs_before <- uniqueN(ASVTAB$ASV)
abundance_before <- sum(ASVTAB$Abundance)

tagjump_events <- sum(ASVTAB$TagJump, na.rm = TRUE)

tagjump_reads <- sum(
  ASVTAB[TagJump == TRUE, Abundance],
  na.rm = TRUE
)


# ------------------------------------------------------------
# Remove tag-jump observations
# ------------------------------------------------------------

filtered_long <- ASVTAB[TagJump == FALSE]

asvs_after <- uniqueN(filtered_long$ASV)
abundance_after <- sum(filtered_long$Abundance)

fwrite(
  filtered_long,
  file.path(out_path, "LIB2_UNCROSS2_filtered_long.tsv"),
  sep = "\t"
)


# ------------------------------------------------------------
# Reconstruct ASV abundance table
# ------------------------------------------------------------

filtered_wide <- dcast(
  filtered_long,
  ASV ~ SampleID,
  value.var = "Abundance",
  fill = 0
)

fwrite(
  filtered_wide,
  file.path(out_path, "LIB2_ASV_table_uncross2.csv")
)


# ------------------------------------------------------------
# Keep matching ASV sequences
# ------------------------------------------------------------

kept_ids <- filtered_wide$ASV

seq_lookup <- data.table(
  ASV = asv_ids,
  Sequence = asv_sequences
)

kept_map <- seq_lookup[match(kept_ids, ASV)]

stopifnot(!any(is.na(kept_map$Sequence)))


# ------------------------------------------------------------
# FASTA with matching ASV IDs
# ------------------------------------------------------------

fasta_lines <- as.vector(
  rbind(
    paste0(">", kept_map$ASV),
    kept_map$Sequence
  )
)

writeLines(
  fasta_lines,
  file.path(out_path, "LIB2_ASVs_uncross2.fasta")
)


# ------------------------------------------------------------
# Save DADA2-style RDS table
# samples x sequences
# ------------------------------------------------------------

sample_cols <- setdiff(names(filtered_wide), "ASV")

mat_ids <- as.matrix(filtered_wide[, ..sample_cols])

rownames(mat_ids) <- filtered_wide$ASV

seq_names <- kept_map$Sequence

seqtab.uncross2 <- t(mat_ids)

colnames(seqtab.uncross2) <- seq_names

saveRDS(
  seqtab.uncross2,
  file.path(out_path, "LIB2_seqtab_uncross2.rds")
)


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summary <- data.table(
  Metric = c(
    "ASVs_before",
    "ASVs_after",
    "ASV_percent_retained",
    "Abundance_before",
    "Abundance_after",
    "Abundance_percent_retained",
    "TagJump_events_removed",
    "TagJump_reads_removed"
  ),

  Value = c(
    asvs_before,
    asvs_after,
    round(100 * asvs_after / asvs_before, 2),
    abundance_before,
    abundance_after,
    round(100 * abundance_after / abundance_before, 2),
    tagjump_events,
    tagjump_reads
  )
)

fwrite(
  summary,
  file.path(out_path, "LIB2_UNCROSS2_summary.tsv"),
  sep = "\t"
)


# ------------------------------------------------------------
# Print results
# ------------------------------------------------------------

cat("\n============================================\n")
cat("LIB2 UNCROSS2 COMPLETED\n")
cat("============================================\n")

cat("Samples:", n_samples, "\n")

cat("\nASVs before:", asvs_before, "\n")
cat("ASVs after:", asvs_after, "\n")

cat(
  "ASV retention:",
  round(100 * asvs_after / asvs_before, 2),
  "%\n"
)

cat("\nAbundance before:", abundance_before, "\n")
cat("Abundance after:", abundance_after, "\n")

cat(
  "Abundance retention:",
  round(100 * abundance_after / abundance_before, 2),
  "%\n"
)

cat("\nTag-jump events removed:", tagjump_events, "\n")
cat("Tag-jump reads removed:", tagjump_reads, "\n")

cat("\nSTOP AND REVIEW BEFORE CONTAMINANT ASSESSMENT.\n")
