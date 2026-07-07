# ============================================================
# Make phylogenetic clusters (alphabetical, nearest-neighbor)
#
# Goal:
#   - Load the working phylogeny (Range ∩ Tree) from:
#       /vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm/raw_data/amphibians/02-phylogeny_pruned_to_working_species_range_x_tree.Rdata
#   - Create clusters of ~30 closest species to a focal species
#   - Process focal species alphabetically (greedy assignment)
#   - Hard constraint: do NOT split "sister species" (tips sharing the same immediate parent)
#
# Outputs:
#   - clusters object (list) + mapping table
#   - CSV: species -> cluster_id
#   - CSV: cluster sizes summary
# ============================================================

suppressPackageStartupMessages({
  library(ape)
})

# -----------------------------
# Paths / parameters
# -----------------------------
tree_rdata_path <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm/raw_data/amphibians/02-phylogeny_pruned_to_working_species_range_x_tree.Rdata"
out_dir         <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm/raw_data/amphibians/clusters"
out_dir_spList <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm/raw_data/amphibians"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

target_cluster_size <- 30L

# -----------------------------
# Load working phylogeny
# -----------------------------
# This file should contain `tree_working` as saved by your harmonization script.
load(tree_rdata_path)

if (!exists("tree_working")) {
  stop("Expected object `tree_working` not found in: ", tree_rdata_path)
}
tree <- tree_working
rm(tree_working)

if (!inherits(tree, "phylo")) stop("Loaded `tree_working` is not a phylo object.")

# Ensure unique tip labels (should already be true after dedup)
if (any(duplicated(tree$tip.label))) {
  dup <- unique(tree$tip.label[duplicated(tree$tip.label)])
  stop("Tree has duplicated tip labels; cannot cluster safely. Example duplicates: ",
       paste(head(dup, 10), collapse = ", "))
}

tips <- sort(tree$tip.label)
Ntip <- length(tips)

cat("Tree tips:", Ntip, "\n")

# -----------------------------
# Sister-pair map (do not split cherries)
# Definition here: two tips that share the same immediate parent node.
# -----------------------------
# Build tip -> parent lookup from edge matrix
# In ape::phylo, tips are 1..Ntip; internal nodes are (Ntip+1)..(Ntip+Nnode)
edge <- tree$edge
parents <- integer(Ntip)
for (t in seq_len(Ntip)) {
  rw <- which(edge[, 2] == t)
  if (length(rw) != 1) stop("Tip ", t, " has ", length(rw), " parent edges; unexpected.")
  parents[t] <- edge[rw, 1]
}

# Identify cherries: internal node with exactly 2 tip children
tip_children_by_parent <- split(seq_len(Ntip), parents)

sister_of <- setNames(rep(NA_character_, Ntip), tree$tip.label)
cherry_parents <- names(tip_children_by_parent)[vapply(tip_children_by_parent, length, integer(1)) == 2L]

for (p in cherry_parents) {
  idx <- tip_children_by_parent[[p]]
  a <- tree$tip.label[idx[1]]
  b <- tree$tip.label[idx[2]]
  sister_of[a] <- b
  sister_of[b] <- a
}

n_cherries <- sum(!is.na(sister_of)) / 2
cat("Identified cherry sister pairs:", n_cherries, "\n")

# Convenience: closure that ensures if a species is included, its cherry-sister is included too (if present)
closure_with_sister <- function(species_vec) {
  species_vec <- unique(species_vec)
  add <- character()
  
  for (s in species_vec) {
    sib <- sister_of[[s]]
    if (!is.na(sib)) add <- c(add, sib)
  }
  unique(c(species_vec, add))
}

# -----------------------------
# Distance matrix (patristic)
# -----------------------------
# cophenetic.phylo returns matrix keyed by tip.label
cat("Computing phylogenetic distances (cophenetic.phylo)...\n")
D <- cophenetic.phylo(tree)

# -----------------------------
# Greedy clustering in alphabetical order
# -----------------------------
assigned <- setNames(rep(FALSE, Ntip), tips)
cluster_id <- setNames(rep(NA_integer_, Ntip), tips)
clusters <- list()

# Helper: add a set of species to a cluster, respecting assignment
add_to_cluster <- function(current, candidates) {
  candidates <- unique(candidates)
  candidates <- candidates[!assigned[candidates]]
  if (length(candidates) == 0) return(current)
  
  # enforce sister closure for any candidate additions
  candidates2 <- closure_with_sister(candidates)
  candidates2 <- unique(candidates2)
  candidates2 <- candidates2[!assigned[candidates2]]
  
  unique(c(current, candidates2))
}

# Main loop: alphabetical focal species
cid <- 0L
for (focal in tips) {
  if (assigned[focal]) next
  
  cid <- cid + 1L
  current <- character()
  
  # Always include focal + its sister (if any)
  current <- add_to_cluster(current, focal)
  
  # Build ordered list of nearest neighbors by distance (excluding self)
  # Use the full D row; then filter to unassigned and not already in cluster
  drow <- D[focal, ]
  drow <- drow[names(drow) != focal]
  ord <- names(sort(drow, decreasing = FALSE))
  
  # Greedily add nearest neighbors until we hit ~target size or exhaust candidates
  i <- 1L
  while (length(current) < target_cluster_size && i <= length(ord)) {
    cand <- ord[i]
    i <- i + 1L
    
    if (assigned[cand]) next
    
    # If candidate is unassigned, propose adding candidate (+ its sister if any),
    # but do not partially add (closure enforced).
    proposed <- add_to_cluster(current, cand)
    
    # If proposed doesn't change, skip
    if (identical(sort(proposed), sort(current))) next
    
    # Accept proposal. Note: this may push cluster size > target due to sister constraint.
    current <- proposed
  }
  
  # If we somehow ended with an empty cluster (should not happen), skip
  if (length(current) == 0) next
  
  # Mark assigned
  assigned[current] <- TRUE
  cluster_id[current] <- cid
  clusters[[as.character(cid)]] <- sort(current)
  
  if (cid %% 25 == 0) cat("Built clusters:", cid, " | assigned:", sum(assigned), "/", Ntip, "\n")
}

cat("\n--- Clustering summary ---\n")
cat("Total clusters:", length(clusters), "\n")
cat("Assigned species:", sum(assigned), " / ", Ntip, "\n", sep = "")

# -----------------------------
# Post checks: ensure no cherry sisters split
# -----------------------------
sister_split <- character()
for (s in names(sister_of)) {
  sib <- sister_of[[s]]
  if (is.na(sib)) next
  if (!is.na(cluster_id[s]) && !is.na(cluster_id[sib]) && cluster_id[s] != cluster_id[sib]) {
    sister_split <- unique(c(sister_split, s, sib))
  }
}
if (length(sister_split) > 0) {
  stop("Sister-splitting detected for: ", paste(head(sort(unique(sister_split)), 20), collapse = ", "),
       "\nThis should not happen; investigate tree edge structure / sister_of mapping.")
} else {
  cat("Check passed: no cherry sister pairs were split.\n")
}

# -----------------------------
# Write outputs
# -----------------------------
mapping_df <- data.frame(
  species = names(cluster_id),
  cluster_id = unname(cluster_id),
  stringsAsFactors = FALSE
)
mapping_df <- mapping_df[order(mapping_df$cluster_id, mapping_df$species), ]

cluster_sizes <- data.frame(
  cluster_id = as.integer(names(clusters)),
  size = vapply(clusters, length, integer(1)),
  stringsAsFactors = FALSE
)
cluster_sizes <- cluster_sizes[order(cluster_sizes$cluster_id), ]

# ============================================================
# Assign human-readable cluster names
# Rule:
#   - Take most frequent genus in cluster
#   - Use first 4 letters of genus
#   - If duplicated across clusters, append numbers (Abav1, Abav2, ...)
# ============================================================

# Helper: extract genus from species name
get_genus <- function(x) sub("_.*$", "", x)

# Step 1: determine base name for each cluster
cluster_base <- sapply(clusters, function(sp) {
  genera <- get_genus(sp)
  tab <- sort(table(genera), decreasing = TRUE)
  
  # Most frequent genus; tie-break alphabetically
  top_freq <- max(tab)
  top_genera <- sort(names(tab)[tab == top_freq])
  
  genus <- top_genera[1]
  base <- substr(genus, 1, 4)
  
  # Title case (Abav instead of abav)
  paste0(toupper(substr(base, 1, 1)), tolower(substr(base, 2, 4)))
})

# Step 2: resolve duplicates by numbering
cluster_name <- character(length(cluster_base))
names(cluster_name) <- names(cluster_base)

used <- list()

for (cid in names(cluster_base)) {
  base <- cluster_base[[cid]]
  
  if (!base %in% names(used)) {
    used[[base]] <- 1
    cluster_name[[cid]] <- paste0(base, "1")
  } else {
    used[[base]] <- used[[base]] + 1
    cluster_name[[cid]] <- paste0(base, used[[base]])
  }
}

# Step 3: attach names to mapping + sizes
cluster_name_df <- data.frame(
  cluster_id   = as.integer(names(cluster_name)),
  cluster_name = unname(cluster_name),
  stringsAsFactors = FALSE
)

# Merge into species → cluster map
mapping_df <- merge(mapping_df, cluster_name_df,
                    by = "cluster_id", all.x = TRUE)
mapping_df <- mapping_df[order(mapping_df$cluster_id, mapping_df$species), ]

# Merge into cluster sizes
cluster_sizes <- merge(cluster_sizes, cluster_name_df,
                       by = "cluster_id", all.x = TRUE)
cluster_sizes <- cluster_sizes[order(cluster_sizes$cluster_id), ]

# ------------------------------------------------------------
# Write updated outputs
# ------------------------------------------------------------
write.csv(mapping_df,
          file = file.path(out_dir, "cluster_map_species_to_cluster_named.csv"),
          row.names = FALSE)

write.csv(cluster_sizes,
          file = file.path(out_dir, "cluster_sizes_named.csv"),
          row.names = FALSE)

# Update saved RData with names
names(clusters) <- cluster_name[match(names(clusters), names(cluster_name))]
save(clusters, 
     file = file.path(out_dir_spList, "spList.Rdata"))
save(mapping_df, cluster_sizes, cluster_name,
     file = file.path(out_dir, "clusters_range_x_tree_nearest30_named.Rdata"))

# Quick sanity print
cat("\n--- Cluster naming preview ---\n")
print(head(cluster_sizes[, c("cluster_id", "cluster_name", "size")], 10))


cat("\nSaved:\n",
    "- ", file.path(out_dir, "cluster_map_species_to_cluster_id.csv"), "\n",
    "- ", file.path(out_dir, "cluster_sizes_summary.csv"), "\n",
    "- ", file.path(out_dir, "clusters_range_x_tree_nearest30_greedy_alphabetical.Rdata"), "\n",
    sep = "")

# Optional: print quick distribution
cat("\nCluster size distribution (min/median/max): ",
    min(cluster_sizes$size), " / ",
    median(cluster_sizes$size), " / ",
    max(cluster_sizes$size), "\n", sep = "")
