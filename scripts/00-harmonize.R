# ============================================================
# Amphibians: harmonize GBIF + expert ranges + phylogeny
# and build a working species list (Range ∩ Tree; GBIF is a flag)
#
# Outputs:
#   - Harmonized + deduplicated phylogeny (RData)
#   - Working species list (Range ∩ Tree) (CSV + TXT)
#   - Pruned phylogeny to working species (RData)
#   - Master species list with presence flags (CSV)
#   - Diagnostics CSVs (missing/ambiguous/coverage)
# ============================================================

suppressPackageStartupMessages({
  library(ape)
})

# -----------------------------
# Paths
# -----------------------------
amph_data_path <- "/vast/palmer/pi/jetz/newt_scratch/standard_occurrences/v1.2/harmonized/amphibians"
amph_tax_path <- "/gpfs/gibbs/pi/jetz/data/species_datasets/mol_taxonomy/MOL_AmphibiaTaxonomy_v2.2_LF.csv"
amph_tree_path <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm/amphibia_mcc_tree_100trees.phy"

# Expert range path: /gpfs/gibbs/pi/jetz/data/species_datasets/rangemaps but using my copy
amph_range_path_local <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm/expert_ranges"

# Where to write checks/outputs (optional checks)
checks_dir <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm/checks"
tree_out_dir <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm/raw_data/amphibians"
dir.create(checks_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tree_out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Read in data
# -----------------------------
amph_tax <- read.csv(amph_tax_path, stringsAsFactors = FALSE)

# Standardize taxonomy to file-name convention (underscores)
amph_tax$Accepted <- gsub(" ", "_", amph_tax$Accepted)
amph_tax$Synonym <- gsub(" ", "_", amph_tax$Synonym)

# Tree
tree_raw <- ape::read.tree(amph_tree_path)

# GBIF species (from harmonized occurrences directory)
gbif_species <- dir(amph_data_path)
gbif_species <- gsub("\\.csv$", "", gbif_species)

# Range species (from YOUR local copy directory)
range_species <- dir(amph_range_path_local)
range_species <- gsub("\\.gpkg$", "", range_species)

# -----------------------------
# Quick checks: coverage vs accepted names
# -----------------------------
gbif_unmapped <- setdiff(gbif_species, amph_tax$Accepted) # ideally empty if upstream harmonized
range_unmapped <- setdiff(range_species, amph_tax$Accepted)

cat("GBIF unmapped vs Accepted:", length(gbif_unmapped), "\n")
cat("Range unmapped vs Accepted:", length(range_unmapped), "\n")

write.csv(data.frame(gbif_unmapped = sort(unique(gbif_unmapped))),
  file = file.path(checks_dir, "00-gbif_unmapped_vs_accepted.csv"),
  row.names = FALSE
)

write.csv(data.frame(range_unmapped = sort(unique(range_unmapped))),
  file = file.path(checks_dir, "00-range_unmapped_vs_accepted.csv"),
  row.names = FALSE
)

# ============================================================
# 0) Rename expert range files using Synonym -> Accepted mapping
# ============================================================

mark_df <- amph_tax[, c("Accepted", "Synonym")]
mark_df$renamed <- FALSE

missing <- character() # names not in Accepted or Synonym
ambiguous <- character() # synonym maps to >1 accepted OR destination already exists OR mv failure
renamed_n <- 0L
accepted_n <- 0L
skipped_exists_n <- 0L
missing_file_n <- 0L

for (i in seq_along(range_species)) {
  range_fsps <- range_species[i]

  # If already accepted, nothing to do
  if (range_fsps %in% mark_df$Accepted) {
    accepted_n <- accepted_n + 1L
    next
  }

  # Find synonym row(s)
  rw <- which(range_fsps == mark_df$Synonym)

  # Not in synonyms either -> truly missing from taxonomy mapping
  if (length(rw) == 0) {
    missing <- c(missing, range_fsps)
    next
  }

  # Multiple matches -> ambiguous synonym
  if (length(unique(mark_df[rw, "Accepted"])) > 1) {
    ambiguous <- c(ambiguous, range_fsps)
    next
  }

  # Exactly the same match -> rename file
  new_name <- unique(mark_df$Accepted[rw])

  src <- file.path(amph_range_path_local, paste0(range_fsps, ".gpkg"))
  dst <- file.path(amph_range_path_local, paste0(new_name, ".gpkg"))

  # If source file doesn't exist in your local copy, record it
  if (!file.exists(src)) {
    missing_file_n <- missing_file_n + 1L
    missing <- c(missing, range_fsps)
    next
  }

  # Avoid overwriting an existing destination
  if (file.exists(dst)) {
    skipped_exists_n <- skipped_exists_n + 1L
    ambiguous <- c(ambiguous, range_fsps)
    next
  }

  cmd <- sprintf("mv %s %s", shQuote(src), shQuote(dst))
  status <- system(cmd)

  if (identical(status, 0L)) {
    mark_df$renamed[rw] <- TRUE
    renamed_n <- renamed_n + 1L
  } else {
    # mv failed for some reason (permissions, etc.)
    ambiguous <- c(ambiguous, range_fsps)
  }
}

cat("\n--- Range rename summary ---\n")
cat("Already Accepted (no change):", accepted_n, "\n")
cat("Renamed via Synonym -> Accepted:", renamed_n, "\n")
cat("Missing (not in Accepted/Synonym):", length(unique(missing)), "\n")
cat("Ambiguous / skipped:", length(unique(ambiguous)), "\n")
cat("Skipped because destination exists:", skipped_exists_n, "\n")
cat("Missing source files in local dir:", missing_file_n, "\n")

write.csv(data.frame(missing_name = sort(unique(missing))),
  file = file.path(checks_dir, "00-harmonize_missing_expert_range_species.csv"),
  row.names = FALSE
)

write.csv(data.frame(ambiguous_name = sort(unique(ambiguous))),
  file = file.path(checks_dir, "00-harmonize_ambiguous_range_species.csv"),
  row.names = FALSE
)

# Verify: what remains in your local directory that is still not in Accepted?
range_species_after <- gsub("\\.gpkg$", "", dir(amph_range_path_local))
still_unmapped_ranges <- setdiff(range_species_after, amph_tax$Accepted)

cat("\nRemaining range filenames not in Accepted after renaming:", length(still_unmapped_ranges), "\n")
if (length(still_unmapped_ranges) > 0) {
  write.csv(data.frame(still_unmapped = sort(still_unmapped_ranges)),
    file = file.path(checks_dir, "00-harmonize_still_unmapped_after_renaming.csv"),
    row.names = FALSE
  )
  print(head(still_unmapped_ranges, 25))
}

# ============================================================
# 1) Harmonize phylogeny tip labels to Accepted taxonomy
#    (translate tip labels in-memory; no file renames)
# ============================================================

# Helper: standardize common tip label formatting
standardize_tip <- function(x) {
  x <- gsub(" ", "_", x)
  x <- gsub("-", "_", x)
  x <- gsub("\\.", "_", x)
  x
}

tree <- tree_raw
tree$tip.label <- standardize_tip(tree$tip.label)

accepted_set <- unique(amph_tax$Accepted)

# Synonym -> Accepted mapping (may be 1-to-many)
syn2acc_list <- split(amph_tax$Accepted, amph_tax$Synonym)

tip0 <- tree$tip.label
tip_new <- tip0

tree_tip_missing <- character() # not in Accepted and not in Synonym
tree_tip_ambiguous <- character() # synonym maps to >1 accepted
tree_tip_changed <- character() # tips that were translated
tree_tip_nochange <- character() # tips already accepted

for (i in seq_along(tip0)) {
  t <- tip0[i]

  if (t %in% accepted_set) {
    tree_tip_nochange <- c(tree_tip_nochange, t)
    next
  }

  if (t %in% names(syn2acc_list)) {
    candidates <- unique(syn2acc_list[[t]])
    if (length(candidates) == 1) {
      tip_new[i] <- candidates[1]
      tree_tip_changed <- c(tree_tip_changed, t)
    } else {
      tree_tip_ambiguous <- c(tree_tip_ambiguous, t)
      # keep original label to avoid incorrect mapping
    }
  } else {
    tree_tip_missing <- c(tree_tip_missing, t)
  }
}

tree$tip.label <- tip_new

# Duplicates created by translation (two tips mapping to same Accepted)
dup_after <- tree$tip.label[duplicated(tree$tip.label)]
dup_after <- sort(unique(dup_after))

cat("\n--- Phylogeny harmonization summary ---\n")
cat("Tips already Accepted:", length(unique(tree_tip_nochange)), "\n")
cat("Tips translated via Synonym->Accepted:", length(unique(tree_tip_changed)), "\n")
cat("Tips missing from taxonomy (not Accepted/Synonym):", length(unique(tree_tip_missing)), "\n")
cat("Tips ambiguous (synonym->multiple Accepted):", length(unique(tree_tip_ambiguous)), "\n")
cat("Duplicate Accepted labels created after translation:", length(dup_after), "\n")

write.csv(data.frame(tree_tip_missing = sort(unique(tree_tip_missing))),
  file = file.path(checks_dir, "01-phylo_tiplabels_missing_from_taxonomy.csv"),
  row.names = FALSE
)

write.csv(data.frame(tree_tip_ambiguous = sort(unique(tree_tip_ambiguous))),
  file = file.path(checks_dir, "01-phylo_tiplabels_ambiguous_synonyms.csv"),
  row.names = FALSE
)

write.csv(data.frame(duplicate_accepted_labels = dup_after),
  file = file.path(checks_dir, "01-phylo_tiplabels_duplicates_after_translation.csv"),
  row.names = FALSE
)

# Resolve duplicates deterministically by dropping extras (keep first occurrence)
if (length(dup_after) > 0) {
  keep_idx <- !duplicated(tree$tip.label)
  tree_dedup <- drop.tip(tree, tree$tip.label[!keep_idx])
} else {
  tree_dedup <- tree
}

cat("Tree tips before/after dedup:", length(tree$tip.label), "->", length(tree_dedup$tip.label), "\n")

# Save harmonized, deduplicated tree
save(tree_dedup, file = file.path(tree_out_dir, "01-phylogeny_harmonized_dedup.Rdata"))

# ============================================================
# 2) Working species list = Range ∩ Phylogeny
#    (GBIF is tracked, but not required)
# ============================================================

range_species_final <- gsub("\\.gpkg$", "", dir(amph_range_path_local))
tree_species_final <- unique(tree_dedup$tip.label)

working_species <- intersect(range_species_final, tree_species_final)
working_species <- sort(unique(working_species))

cat("\n--- Working species list summary (Range ∩ Tree) ---\n")
cat("GBIF species:", length(unique(gbif_species)), "\n")
cat("Range species (after renaming):", length(unique(range_species_final)), "\n")
cat("Tree species (harmonized, dedup):", length(unique(tree_species_final)), "\n")
cat("Working species (Range ∩ Tree):", length(working_species), "\n")
cat("Working species WITH GBIF data:", sum(working_species %in% gbif_species), "\n")
cat("Working species MISSING GBIF data:", sum(!(working_species %in% gbif_species)), "\n")

# Diagnostics under Range ∩ Tree definition
range_not_in_working <- setdiff(range_species_final, working_species)
tree_not_in_working <- setdiff(tree_species_final, working_species)
gbif_missing_within_working <- setdiff(working_species, gbif_species)

write.csv(data.frame(working_species = working_species),
  file = file.path(checks_dir, "02-working_species_list_range_x_tree.csv"),
  row.names = FALSE
)

writeLines(working_species, con = file.path(checks_dir, "02-working_species_list_range_x_tree.txt"))

write.csv(data.frame(range_not_in_working = sort(unique(range_not_in_working))),
  file = file.path(checks_dir, "02-range_not_in_working_range_x_tree.csv"),
  row.names = FALSE
)

write.csv(data.frame(tree_not_in_working = sort(unique(tree_not_in_working))),
  file = file.path(checks_dir, "02-tree_not_in_working_range_x_tree.csv"),
  row.names = FALSE
)

write.csv(data.frame(gbif_missing_within_working = sort(unique(gbif_missing_within_working))),
  file = file.path(checks_dir, "02-gbif_missing_within_working_range_x_tree.csv"),
  row.names = FALSE
)

# Prune tree to working species (Range ∩ Tree)
tree_working <- drop.tip(tree_dedup, setdiff(tree_dedup$tip.label, working_species))
save(tree_working, file = file.path(tree_out_dir, "02-phylogeny_pruned_to_working_species_range_x_tree.Rdata"))

# ============================================================
# 3) Master species list with presence flags:
#    columns: expert_range, gbif, phylogeny, in_working
#    Master namespace = union(Range, Tree, GBIF)
# ============================================================

master_species <- sort(unique(c(range_species_final, tree_species_final, gbif_species)))

master_df <- data.frame(
  species = master_species,
  expert_range = master_species %in% range_species_final,
  gbif = master_species %in% gbif_species,
  phylogeny = master_species %in% tree_species_final,
  stringsAsFactors = FALSE
)

# Working definition: Range ∩ Tree (GBIF not required)
master_df$in_working <- master_df$expert_range & master_df$phylogeny

# Useful summary counts
master_df$missing_count_all3 <- 3L - (master_df$expert_range + master_df$gbif + master_df$phylogeny)
master_df$missing_count_working_def <- 2L - (master_df$expert_range + master_df$phylogeny)

# Write master list
write.csv(master_df,
  file = file.path(checks_dir, "03-master_species_list_presence_flags.csv"),
  row.names = FALSE
)

# Additional diagnostics from master list
only_tree <- master_df$species[master_df$phylogeny & !master_df$expert_range]
only_range <- master_df$species[master_df$expert_range & !master_df$phylogeny]
in_range_tree_no_gbif <- master_df$species[master_df$in_working & !master_df$gbif]

write.csv(data.frame(species = sort(unique(only_tree))),
  file = file.path(checks_dir, "03-only_phylogeny_not_range.csv"),
  row.names = FALSE
)

write.csv(data.frame(species = sort(unique(only_range))),
  file = file.path(checks_dir, "03-only_range_not_phylogeny.csv"),
  row.names = FALSE
)

write.csv(data.frame(species = sort(unique(in_range_tree_no_gbif))),
  file = file.path(checks_dir, "03-in_working_missing_gbif.csv"),
  row.names = FALSE
)

cat("\n--- Master list summary ---\n")
cat("Master species (union of Range, Tree, GBIF):", nrow(master_df), "\n")
cat("In phylogeny:", sum(master_df$phylogeny), "\n")
cat("In expert ranges:", sum(master_df$expert_range), "\n")
cat("In GBIF:", sum(master_df$gbif), "\n")
cat("Working (Range ∩ Tree):", sum(master_df$in_working), "\n")
cat("Working but missing GBIF:", sum(master_df$in_working & !master_df$gbif), "\n")

cat("\nSaved key outputs:\n",
  "- Harmonized dedup tree:  ", file.path(tree_out_dir, "01-phylogeny_harmonized_dedup.Rdata"), "\n",
  "- Pruned working tree:    ", file.path(tree_out_dir, "02-phylogeny_pruned_to_working_species_range_x_tree.Rdata"), "\n",
  "- Working species list:   ", file.path(checks_dir, "02-working_species_list_range_x_tree.csv"), "\n",
  "- Master species flags:   ", file.path(checks_dir, "03-master_species_list_presence_flags.csv"), "\n",
  "- Diagnostics in:         ", checks_dir, "\n",
  sep = ""
)
