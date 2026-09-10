library(dada2)

base <- "/home/sofija/INTERNSHIP/samples/Capreolus_capreolus/Plants_ITS"

raw_path <- file.path(base, "dada2_out", "raw")
nochim_path <- file.path(base, "dada2_out", "nochim")

dir.create(nochim_path, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load raw ASV table and ASV mapping
# ------------------------------------------------------------

seqtab <- readRDS(
  file.path(raw_path, "LIB1_seqtab_raw.rds")
)

raw_map <- read.csv(
  file.path(raw_path, "LIB1_ASV_sequences_raw.csv"),
  stringsAsFactors = FALSE
)

stopifnot(identical(colnames(seqtab), raw_map$Sequence))

cat("ASVs before chimera removal:", ncol(seqtab), "\n")
cat("Abundance before chimera removal:", sum(seqtab), "\n")


# ------------------------------------------------------------
# Chimera removal
# ------------------------------------------------------------

seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

saveRDS(
  seqtab.nochim,
  file.path(nochim_path, "LIB1_seqtab_nochim.rds")
)


# ------------------------------------------------------------
# Preserve original ASV identifiers
# ------------------------------------------------------------

keep <- match(colnames(seqtab.nochim), raw_map$Sequence)

stopifnot(!any(is.na(keep)))

nochim_ids <- raw_map$ASV_ID[keep]
nochim_sequences <- colnames(seqtab.nochim)

nochim_map <- data.frame(
  ASV_ID = nochim_ids,
  Sequence = nochim_sequences
)

write.csv(
  nochim_map,
  file.path(nochim_path, "LIB1_ASV_sequences_nochim.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# Chimera-filtered ASV abundance table
# ------------------------------------------------------------

nochim_table <- as.data.frame(t(seqtab.nochim))
rownames(nochim_table) <- nochim_ids

write.csv(
  nochim_table,
  file.path(nochim_path, "LIB1_ASV_table_nochim.csv"),
  row.names = TRUE
)


# ------------------------------------------------------------
# FASTA
# ------------------------------------------------------------

fasta_lines <- as.vector(
  rbind(
    paste0(">", nochim_ids),
    nochim_sequences
  )
)

writeLines(
  fasta_lines,
  file.path(nochim_path, "LIB1_ASVs_nochim.fasta")
)


# ------------------------------------------------------------
# Per-sample read tracking
# ------------------------------------------------------------

merged_reads <- rowSums(seqtab)
nonchim_reads <- rowSums(seqtab.nochim)

tracking <- data.frame(
  Sample = rownames(seqtab),
  Merged = merged_reads,
  Non_chimeric = nonchim_reads,
  Percent_retained = round(
    100 * nonchim_reads / merged_reads,
    2
  )
)

write.csv(
  tracking,
  file.path(nochim_path, "LIB1_chimera_tracking.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

raw_asvs <- ncol(seqtab)
nochim_asvs <- ncol(seqtab.nochim)

raw_abundance <- sum(seqtab)
nochim_abundance <- sum(seqtab.nochim)

cat("\n============================================\n")
cat("LIB1 CHIMERA REMOVAL COMPLETED\n")
cat("============================================\n")

cat("ASVs before:", raw_asvs, "\n")
cat("ASVs after:", nochim_asvs, "\n")
cat(
  "ASV retention:",
  round(100 * nochim_asvs / raw_asvs, 2),
  "%\n"
)

cat("\nAbundance before:", raw_abundance, "\n")
cat("Abundance after:", nochim_abundance, "\n")
cat(
  "Abundance retention:",
  round(100 * nochim_abundance / raw_abundance, 2),
  "%\n"
)

cat("\nPer-sample tracking:\n")
print(tracking)

cat("\nSTOP AND REVIEW BEFORE UNCROSS2.\n")

