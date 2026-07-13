# phyloSDM_MOL

Phylogenetic species distribution models (SDMs) for species, built on occurrence data, expert range maps, and a dated phylogeny from [Map of Life](https://mol.org) (MOL). The core model is a **Log-Gaussian Cox Process (LGCP)** fit in Stan (with an optional TMB implementation), where each species' response-to-environment coefficients are given a **phylogenetic Gaussian Process prior** — closely related species are expected to respond similarly to climate/habitat covariates, and the model borrows statistical strength across the tree accordingly.

This pipeline was developed for and run on Yale's McCleary HPC cluster (SLURM + [dSQ](https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/dsq/)), but every step can also be run interactively on a single machine (see [Running locally vs. on an HPC cluster](#running-locally-vs-on-an-hpc-cluster)).

---

## Table of contents

1. [What this project does](#what-this-project-does)
2. [Repository layout](#repository-layout)
3. [Getting the data](#getting-the-data)
4. [Software prerequisites](#software-prerequisites)
5. [Setting up the conda environment](#setting-up-the-conda-environment)
6. [Before you run anything: required edits](#before-you-run-anything-required-edits)
7. [Configuring an experiment (`01-config.yaml`)](#configuring-an-experiment-01-configyaml)
8. [Pipeline walkthrough](#pipeline-walkthrough)
9. [Running locally vs. on an HPC cluster](#running-locally-vs-on-an-hpc-cluster)
10. [Quickstart: running one cluster end-to-end](#quickstart-running-one-cluster-end-to-end)
11. [Output directory reference](#output-directory-reference)
12. [Troubleshooting](#troubleshooting)
13. [Open questions for the project author](#open-questions-for-the-project-author)
14. [Citation](#citation)

---

## What this project does

Given:
- Point occurrence records for a set of species (e.g. from GBIF, harmonized against MOL taxonomy),
- Expert-drawn range maps for those species,
- A time-calibrated phylogeny spanning them, and
- A stack of environmental rasters (temperature, precipitation, cloud cover, vegetation index, topographic ruggedness, elevation — from [CHELSA](https://chelsa-climate.org/) and related products),

the pipeline:

1. Harmonizes taxonomy across all three data sources and builds a "working species list" (species present in **both** the range maps and the phylogeny).
2. Splits that species list into phylogenetically coherent **clusters** of ~30 species each (so each model fit involves a tractable covariance matrix — see [why clustering?](#why-cluster-species)).
3. For each cluster: builds presence/background training data, fits the phylogenetic LGCP model in Stan (or TMB), evaluates predictive performance on held-out data, and produces continuous and thresholded binary range-map predictions per species.

The end product is, per species, a continuous "relative probability of occurrence" raster and a thresholded binary presence/absence raster — i.e. a model-based species range map that leverages phylogenetic relatedness to improve predictions, especially for data-poor species.

---

## Repository layout

This repo's `.gitignore` is deliberately aggressive: **only `scripts/` is tracked in git.** Everything else (raw data, intermediate outputs, model fits, predictions) is generated locally or distributed separately via Zenodo, because most of it is either large binary geospatial data or reproducible from the scripts.

```
phyloSDM_MOL/
├── scripts/            # ✅ tracked in git — all R code, the pipeline config, and conda env files
├── raw_data/            # ❌ not tracked — phylogeny, taxonomy, per-species occurrence CSVs, per-cluster intermediates
│   └── amphibians/
├── expert_ranges/        # ❌ not tracked — per-species expert range maps (.gpkg) + a combined shapefile
├── analysis/            # ❌ not tracked — soft-clip rasters, conditional predictions, evaluation, spatial predictions, binary maps
│   └── soft_clips/
├── data/               # ❌ not tracked — currently empty; not written to by any script (see open questions)
├── res/                # ❌ not tracked — model fit objects (.Rdata) land here
├── jobs/               # ❌ not tracked — SLURM/dSQ job scripts and generated job-array task lists
├── log/                # ❌ not tracked — SLURM job logs, job manifests
└── README.md
```

**Practical implication:** cloning this repo from GitHub gets you *code only*. You must separately obtain `raw_data/`, `expert_ranges/`, and (if you want to skip re-running early pipeline stages) `analysis/` from the Zenodo archive below.

---

## Getting the data

The data underlying this project (harmonized phylogeny, expert range maps, occurrence CSVs, and — depending on what's included — intermediate/output files) is archived on Zenodo:

> **Zenodo DOI:** 10.5281/zenodo.21243412
>
> ```bash
> # once you have the DOI, something like:
> curl -L -o phyloSDM_MOL_data.zip "https://zenodo.org/record/XXXXXXX/files/phyloSDM_MOL_data.zip"
> unzip phyloSDM_MOL_data.zip -d .
> ```

After extraction, you should end up with `raw_data/`, `expert_ranges/`, and (if included) `analysis/` sitting directly alongside `scripts/` at the repo root, matching the [layout above](#repository-layout).

### A separate, large environmental-raster dependency

Independent of the Zenodo archive, the data-generation script (`03-gen_data.R`) expects a directory of global environmental rasters **outside the repo**, at `~/env` (when run locally) — roughly 30 GB, containing:

| File | Layer |
|---|---|
| `CHELSA_bio_1.tif` | Mean annual temperature |
| `CHELSA_bio_4.tif` | Temperature seasonality |
| `CHELSA_bio_13.tif` | Precipitation of wettest month |
| `CHELSA_bio_15.tif` | Precipitation seasonality |
| `cloudCover.tif` | Cloud cover |
| `Annual_EVI.tif` | Annual enhanced vegetation index |
| `TRI.tif` | Topographic ruggedness index |
| `elevation_1KMmean_SRTM.tif` | Elevation |

⚠️ Gather env rasters independently from source

---

## Software prerequisites

- **[Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/)** — used to build the R + geospatial environment.
- **Git**.
- **A C/C++ compiler** — only required if you plan to fit models with `MODEL_TYPE = "TMB"` (see [pipeline step 05](#pipeline-walkthrough)).
- **`cmdstanr` + [CmdStan](https://mc-stan.org/docs/cmdstan-guide/) — required for `MODEL_TYPE = "STAN"` when running locally (off the HPC).** `rstan` (in the conda env) cannot compile Stan models on current macOS: `StanHeaders`/`RcppParallel`'s bundled TBB headers are incompatible with recent Apple Clang (see [Troubleshooting](#troubleshooting)). `05-run_model.R` uses `cmdstanr` instead when running locally, which compiles via CmdStan's own build system and avoids this entirely. Install with `install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))` then `cmdstanr::install_cmdstan()`, from whichever R installation you run the pipeline with locally (**not necessarily the conda env** — see below). On the HPC, `05-run_model.R` keeps using `rstan` (unaffected by this issue), so nothing extra is needed there.
- **(HPC only)** SLURM (`sbatch`) and Yale's [dSQ](https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/dsq/) module — used for submitting per-cluster jobs as array jobs. Not needed to run the pipeline on a single machine.

---

## Setting up the conda environment

The environment specification lives in `scripts/`. There are two variants:

- **`scripts/phylo-sdms_wf.yaml`** — the canonical environment, exported from the Linux HPC system this project runs on (`linux-64`, exact build-hash pins). Use this if you're on Linux or on the HPC cluster itself.
- **`scripts/phylo-sdms_wf.osx.v3.yaml`** — a macOS (Apple Silicon / `osx-arm64`) build, re-solved from the same package *names* since the Linux build hashes don't exist for macOS. Package versions may differ slightly from the Linux environment as a result (see [Troubleshooting](#troubleshooting) for why).

```bash
# Linux / HPC:
conda env create -f scripts/phylo-sdms_wf.yaml

# macOS (Apple Silicon):
conda env create -f scripts/phylo-sdms_wf.osx.v3.yaml

# Activate (both variants use the same env name):
conda activate phylo-sdms_wf
```

This environment includes R 4.3, `rstan`/`StanHeaders`, `terra`, `sf`, `ape`, GDAL/GEOS/PROJ, and the rest of the R package stack the pipeline depends on. **`TMB` is not included** in either environment file — install it separately (`install.packages("TMB")` inside R, from the activated env) only if you intend to use `MODEL_TYPE = "TMB"`. **`cmdstanr`/`posterior` are also not included** — needed for `MODEL_TYPE = "STAN"` locally on macOS (see [Software prerequisites](#software-prerequisites)); install them the same way, or run the STAN path from a separate R installation that already has `cmdstanr` set up.

---

## Before you run anything: required edits

This codebase was written for one specific HPC user/path layout and has **not yet been generalized** for other users or machines. Before running locally, you will need to make the following edits:

### 1. The `HPC` environment variable defaults the *wrong way*

Every numbered script (`03` through `15`) branches on:

```r
HPC <- Sys.getenv("HPC")
if (HPC != "FALSE") {
  # ... HPC paths ...
} else {
  # ... local paths ...
}
```

`Sys.getenv("HPC")` returns `""` (empty string) if the variable is unset — and `"" != "FALSE"` is `TRUE`. **This means the scripts take the HPC code path by default unless you explicitly set `HPC=FALSE`.** Always run local jobs like this:

```bash
HPC=FALSE Rscript scripts/03-gen_data.R ...
```

### 2. Hardcoded root paths

On the "local" branch, every script sets:

```r
root <- "~/phyloSDM_MOL"
```

This assumes you've cloned the repo to exactly `~/phyloSDM_MOL`. If you cloned it somewhere else, edit the `root <- "~/phyloSDM_MOL"` line in each of: `03-gen_data.R`, `03-soft_clips.R`, `04-data_indices.R`, `04-pkg_data.R`, `05-run_model.R`, `10-cond_pred.R`, `11-pkg_test.R`, `12-eval.R`, `13-continuous_pred.R`, `14-thresholding.R`, `15-binary_prediction.R`, `02-build_jobs.R`.

⚠️ **Gotcha:** if a directory with an old/different hardcoded name (e.g. a stale prior checkout) already exists on your machine at that path, R will happily read/write there instead — silently, with no error, since the path *does* exist. Symptoms look like missing data (e.g. "soft clip files found for 0/N observations") rather than a file-not-found crash. If output looks wrong, double check `root` actually resolves to the repo you think it does.

On the "HPC" branch, `root` and the environmental-raster path `epath` are hardcoded to the author's Yale scratch/PI space (`/vast/palmer/pi/jetz/ss4224/...`) — irrelevant off the Yale cluster, but relevant if you're a lab member reusing this on McCleary under a different account.

### 3. `00-harmonize.R` and `00-makeClusters.R` have no local branch at all

Unlike every other script, these two have **only** hardcoded HPC-style absolute paths (tree file, MOL taxonomy CSV, expert range directory, output directories). You must manually edit the path variables at the top of each script to point at your own copies before running them locally. You only need to run these two scripts if you're building a working species list / cluster assignment from scratch (e.g. a different clade or an updated taxonomy) — if you're reproducing the published amphibian analysis, the outputs of these two scripts (`raw_data/amphibians/spList.Rdata` etc.) should already be included in the Zenodo data.

### 4. `01-config.yaml`'s default `script_name` doesn't exist in this repo

```yaml
script:
  script_name: 03-soft_clips.R
```


### 5. Personal identifiers

`01-config.yaml` (`MAIL_USER`) and `jobs/flex_job.sh` (`--mail-user`) are set to the original author's email — update if you're submitting your own SLURM jobs.

---

## Configuring an experiment (`01-config.yaml`)

This file drives `02-build_jobs.R`, which expands it into a grid of jobs (one per cluster × repetition, etc.):

```yaml
data:
  exp_root: v0           # top-level experiment name — namespaces output files
  exp_id: sub1000         # sub-experiment / run ID
  cluster_name: ALL        # a specific cluster name, or ALL to run every cluster
  raw_data: amphibians      # dataset name — matches raw_data/<raw_data>/
  repno: 1              # which train/test repetition, or ALL
  focal_sp: ALL           # ALL (all species as one cluster job), or run_all_sep (one job per species)
  autocommit: FALSE        # auto git-commit outputs (leave FALSE unless you know you want this)
  nrep: 1               # total number of train/test repetitions

script:
  script_name: 03-gen_data.R    # which pipeline script this config is generating jobs for

model_specs:
  model_type: "STAN"        # "STAN" or "TMB"
  model_name: ["LGCP_background"]  # model variant(s) to run

hpc_specs:                    # SLURM/dSQ resource requests — HPC only
  dsq: TRUE
  MEM: '50G'
  TIME: '10:00:00'
  ...
```

Every numbered script from `03` onward expects **the same 9 positional command-line arguments**, in this order:

```
EXP_ROOT  EXP_ID  DATASET  CLUSTER  FSP  REPNO  NREP  MODEL_TYPE  MODEL_NAME
```

`02-build_jobs.R` generates these calls for you (as a dSQ task list or individual `sbatch` scripts) from `01-config.yaml`; you can also invoke any script directly with these 9 arguments for testing (see [Quickstart](#quickstart-running-one-cluster-end-to-end)).

---

## Pipeline walkthrough

| # | Script | Purpose | Reads | Writes | How it's run |
|---|--------|---------|-------|--------|---------------|
| 00 | `00-harmonize.R` | Harmonize GBIF occurrences, expert-range filenames, and phylogeny tip labels to one taxonomy; build the "working species list" (Range ∩ Tree) | MOL taxonomy CSV, MCC tree (`.phy`), expert range dir, harmonized GBIF occurrences | `01-phylogeny_harmonized_dedup.Rdata`, `02-phylogeny_pruned_to_working_species_range_x_tree.Rdata`, diagnostic CSVs | `Rscript 00-harmonize.R` (no CLI args; hardcoded paths — see [required edits](#3-00-harmonizer-and-00-makeclustersr-have-no-local-branch-at-all)) |
| 00 | `00-makeClusters.R` | Greedily partition the phylogeny into ~30-species clusters, never splitting cherry sister-pairs; names clusters after their dominant genus (e.g. `Chal1`, `Kass1`, `Rani1`) | pruned phylogeny `.Rdata` | `spList.Rdata` (cluster → species list), cluster map/size CSVs | `Rscript 00-makeClusters.R` (no CLI args) |
| 02 | `02-build_jobs.R` | Expand `01-config.yaml` (clusters × reps × focal species) into a job grid; write a dSQ task list or per-job `sbatch` scripts | `01-config.yaml`, `spList.Rdata` | `jobs/job_array_*.txt` or per-cluster `.sh` files, job manifest in `log/` | `Rscript 02-build_jobs.R` |
| 03 | `03-gen_data.R` | Per cluster: rasterize presences, generate background/pseudo-absence points, extract & scale environmental covariates | species CSVs, phylogeny, expert ranges, CHELSA rasters, land mask (`world.Rdata`) | `raw_data/<dataset>/<cluster>/<cluster>_run_files.Rdata`, `_extent.csv`, `_env_scales.csv`, `_missing_species.csv` | SLURM/dSQ array, one task per cluster (`FSP="ALL"`) |
| 03 | `03-soft_clips.R` | Per cluster: build a per-species "soft clip" raster — a logistic distance-decay mask around each species' expert range, used later as a model log-offset | `spList.Rdata`, `<cluster>_extent.csv` (from `03-gen_data.R`), expert ranges, CHELSA template | `analysis/soft_clips/<species>_soft_clip.tif` | SLURM/dSQ array, one task per cluster — **must run after `03-gen_data.R`** |
| 04 | `04-data_indices.R` | Build 70/30 train/test row-index splits (× `nrep` repetitions) for every cluster | `<cluster>_run_files.Rdata` (all clusters) | `<cluster>_indices.Rdata` | Not designed for dSQ — run as a single interactive/batch job covering all clusters |
| 04 | `04-pkg_data.R` | Package the training split into a Stan/TMB-ready `model_data` list: phylogenetic distance matrix, X/y matrices, quadratic terms, soft-clip log-offset | `_run_files.Rdata`, `_indices.Rdata`, soft clips | `<exp_root>_<exp_id>_<cluster>_<fsp>_rep_<repno>_model_data.Rdata` | SLURM/dSQ array, one task per cluster × rep |
| 05 | `05-run_model.R` | Fit the LGCP model. `MODEL_TYPE="STAN"`: on HPC, `rstan::stan()` (model from `05-stan_model.R`); running locally, `cmdstanr` (model from `LGCP_background.stan`) instead, since `rstan` can't compile on current macOS (see [Troubleshooting](#troubleshooting)) — the branch is automatic, based on the same `HPC` env var every script already checks. `MODEL_TYPE="TMB"`: compiles `05_lgcp_corrected.cpp` on either platform. | `model_data.Rdata` | `res/<exp_root>_<exp_id>_<cluster>_<fsp>_rep_<repno>_<model_type>_fit.Rdata` — for local STAN fits, `result$posterior` holds plain-array posterior draws (`B`, `alpha`, `rho`, `sigma_f`) plus `result$diagnostics` (divergences, Rhat, ESS), since a live `cmdstanr` fit object can't be reloaded in a later session; HPC STAN fits still save the full `result$fit` `stanfit` object as before | SLURM/dSQ array, one task per cluster × rep |
| 05 | `05-stan_model.R` | Defines the Stan program string for `LGCP_background`, used only on the **HPC/rstan** path: Poisson-log-link LGCP with a phylogenetic-GP prior (squared-exponential kernel on patristic distance) on species coefficients | — | — | `source()`d by `05-run_model.R`'s HPC branch; rstan JIT-compiles it, no manual `stanc` step |
| 05 | `LGCP_background.stan` | The same model as `05-stan_model.R`, as a real `.stan` file for the **local/cmdstanr** path (current Stan array syntax; `offset` renamed `log_offset` since it's now a reserved word) | — | Compiled binary cached next to the source by `cmdstanr::cmdstan_model()` (same caching pattern as TMB's `.o`/`.so`, skipped if source unchanged) | Compiled by `05-run_model.R`'s local branch |
| 05 | `05_lgcp_corrected.cpp` | TMB (C++ autodiff) reimplementation of the same model, used only when `MODEL_TYPE="TMB"` | — | Compiled `.o`/`.so`, gitignored (platform-specific, rebuilt via `TMB::compile()` — see [Troubleshooting](#troubleshooting)) | Compiled via `TMB::compile()` inside `05-run_model.R` if missing/stale |
| 10 | `10-cond_pred.R` | Per species: leave-one-out phylogenetic-GP **conditional** posterior prediction of its coefficients, using the fitted coefficients of the *other* species in (and beyond) its cluster | fit `.Rdata`, training `model_data.Rdata` | `analysis/<exp_root>/cond_pred/*.Rdata` | SLURM/dSQ array or interactive |
| 11 | `11-pkg_test.R` | Package the held-out **test** split, mirroring `04-pkg_data.R` | `_run_files.Rdata`, `_indices.Rdata`, training `model_data.Rdata` | `*_test_data.Rdata` | SLURM/dSQ array |
| 12 | `12-eval.R` | Evaluate both the raw fit and the conditional prediction against test data (AUC, sensitivity/specificity, confusion matrix via `pROC`); select the best model per species | `test_data.Rdata`, fit, `cond_pred.Rdata` | `analysis/<exp_root>/eval/*_eval_{model,cond,combined}.Rdata`, `*_best_model.csv`, per-species prediction CSVs | SLURM/dSQ array — **needs outputs of both `10` and `11`** |
| 13 | `13-continuous_pred.R` | Generate a continuous "relative probability of occurrence" raster per species, using whichever model (`10` vs. `05`) `12-eval.R` selected, plus the soft-clip offset | fit, `best_model.csv`, CHELSA rasters, soft clips | `analysis/<exp_root>/spatial_pred/*_relprob.tif` | SLURM/dSQ array |
| 14 | `14-thresholding.R` | Find each species' ROC/AUC-optimal (Youden) binarization threshold from test data + the continuous prediction raster | `test_data.Rdata`, `*_relprob.tif` | `analysis/<exp_root>/eval/thresholds/*_thresholds.{csv,Rdata}` | SLURM/dSQ array |
| 15 | `15-binary_prediction.R` | Apply the threshold to produce a binary presence/absence map per species (optionally cropped to the expert range extent) | `*_thresholds.csv`, `*_relprob.tif`, expert ranges | `analysis/<exp_root>/binary_maps/*_binary.tif`, summary CSV | SLURM/dSQ array |

Dependency order in one line: **00 → 02 → 03 (gen_data, then soft_clips) → 04 (indices, then pkg_data) → 05 → {10, 11 in parallel} → 12 → 13 → 14 → 15.**

### Why cluster species?

The phylogenetic-GP prior requires a dense J×J covariance matrix (Cholesky-decomposed once per model fit) over the J species in a fit — an O(J³) operation. Fitting the entire amphibian tree (thousands of tips) in a single Stan/TMB run is computationally infeasible, so `00-makeClusters.R` partitions species into phylogenetically coherent groups of a tractable size (~30). `10-cond_pred.R` is what lets a species "borrow" information from close relatives *outside* its own fitted cluster, via the phylogenetic kernel, at prediction time.

---

## Running locally vs. on an HPC cluster

Every script (except `00-harmonize.R`/`00-makeClusters.R`) checks the `HPC` environment variable. **Remember it defaults to the HPC path unless explicitly overridden** (see [required edits](#1-the-hpc-environment-variable-defaults-the-wrong-way)):

```bash
# Local:
HPC=FALSE Rscript scripts/03-gen_data.R v0 test amphibians Chal1 ALL 1 1 STAN LGCP_background

# On McCleary (or another SLURM+dSQ cluster), via a generated job array:
module load miniconda dSQ
conda activate phylo-sdms_wf
dsq --job-file jobs/job_array_<exp_id>_for_03-gen_data.txt --mem-per-cpu 50G -t 10:00:00 --mail-type ALL --submit
```

`02-build_jobs.R` writes the dSQ task-list file for you from `01-config.yaml`, but does not itself call `dsq` — actually submitting the array job is a manual step on the cluster.

---

## Quickstart: running one cluster end-to-end

This assumes you've already: cloned the repo to `~/phyloSDM_MOL`, extracted the Zenodo data into place, and built and activated the conda environment. If you cloned to a different path, see [required edits, item 2](#2-hardcoded-root-paths) first.

```bash
conda activate phylo-sdms_wf

# 0. (Skip if spList.Rdata etc. are already in the Zenodo data)
Rscript scripts/00-harmonize.R
Rscript scripts/00-makeClusters.R

# 3. Generate training data + soft-clip rasters for one cluster (e.g. "Chal1")
HPC=FALSE Rscript scripts/03-gen_data.R   v0 quickstart amphibians Chal1 ALL 1 1 STAN LGCP_background
HPC=FALSE Rscript scripts/03-soft_clips.R v0 quickstart amphibians Chal1 ALL 1 1 STAN LGCP_background

# 4. Train/test split (all clusters, one pass) + package training data for this cluster
HPC=FALSE Rscript scripts/04-data_indices.R v0 quickstart amphibians ALL   ALL 1 1 STAN LGCP_background
HPC=FALSE Rscript scripts/04-pkg_data.R     v0 quickstart amphibians Chal1 ALL 1 1 STAN LGCP_background

# 5. Fit the model
HPC=FALSE Rscript scripts/05-run_model.R v0 quickstart amphibians Chal1 ALL 1 1 STAN LGCP_background

# 10/11. Conditional prediction + test-data packaging (order doesn't matter between these two)
HPC=FALSE Rscript scripts/10-cond_pred.R v0 quickstart amphibians Chal1 ALL 1 1 STAN LGCP_background
HPC=FALSE Rscript scripts/11-pkg_test.R  v0 quickstart amphibians Chal1 ALL 1 1 STAN LGCP_background

# 12-15. Evaluate, predict continuously, threshold, binarize
HPC=FALSE Rscript scripts/12-eval.R              v0 quickstart amphibians Chal1 ALL 1 1 STAN LGCP_background
HPC=FALSE Rscript scripts/13-continuous_pred.R   v0 quickstart amphibians Chal1 ALL 1 1 STAN LGCP_background
HPC=FALSE Rscript scripts/14-thresholding.R      v0 quickstart amphibians Chal1 ALL 1 1 STAN LGCP_background
HPC=FALSE Rscript scripts/15-binary_prediction.R v0 quickstart amphibians Chal1 ALL 1 1 STAN LGCP_background
```

`Chal1` is a small *Chalcorana*-dominated cluster and a reasonable first cluster to test the pipeline on. Swap in any cluster name from `raw_data/amphibians/clusters/cluster_map_species_to_cluster_named.csv` (or `ALL`, where a script supports it) to run others.

---

## Output directory reference

After a full run, expect:

```
raw_data/amphibians/<cluster>/
  <cluster>_run_files.Rdata          # from 03-gen_data.R
  <cluster>_extent.csv
  <cluster>_env_scales.csv
  <cluster>_indices.Rdata            # from 04-data_indices.R
  <cluster>_rep_<n>_model_data.Rdata # from 04-pkg_data.R

analysis/
  soft_clips/<species>_soft_clip.tif       # from 03-soft_clips.R
  <exp_root>/cond_pred/...                 # from 10-cond_pred.R
  <exp_root>/eval/...                      # from 12-eval.R, 14-thresholding.R
  <exp_root>/spatial_pred/*_relprob.tif    # from 13-continuous_pred.R
  <exp_root>/binary_maps/*_binary.tif      # from 15-binary_prediction.R

res/
  <exp_root>_<exp_id>_<cluster>_<fsp>_rep_<n>_<model_type>_fit.Rdata  # from 05-run_model.R
```

`data/` remains empty throughout — no script in this pipeline currently reads from or writes to it (`⚠️ CONFIRM` — see [Open questions](#open-questions-for-the-project-author)).

---

## Troubleshooting

**`conda env create` fails with `PackagesNotFoundError` on macOS.** `scripts/phylo-sdms_wf.yaml` was exported from the Linux HPC environment with exact build-hash pins (e.g. `gcc_impl_linux-64`, `libgcc-ng`) that don't exist for macOS. Use `scripts/phylo-sdms_wf.osx.v3.yaml` instead — it targets the same package versions where possible but lets conda re-solve build hashes and, in a few cases (the GDAL/GEOS/PROJ/sqlite/AWS-SDK cluster), lets versions float where the exact Linux pins had no macOS equivalent.

**TMB won't compile / `05_lgcp_corrected.so` fails to load, or `make: Nothing to be done` and the model behaves oddly.** `.o`/`.so` files are platform-specific and gitignored (not committed) for exactly this reason — if you have stale ones lying around (e.g. copied from a Linux machine, or built earlier under a different toolchain), delete `scripts/05_lgcp_corrected.o` and `.so` and let `TMB::compile()` (called from `05-run_model.R`) rebuild them fresh. This requires `TMB` to be installed (not in the conda env by default) and a working C++ compiler (Xcode Command Line Tools on macOS). `make` only recompiles when the source is newer than the existing `.o` — a stale binary sitting there is silently treated as "up to date."

**`rstan::stan()` fails to compile with `Error in compileCode(...)` on macOS**, e.g. `a template argument list is expected after a name prefixed by the template keyword` or `'tbb/tbb_stddef.h' file not found` or `in-class initializer for static data member is not a constant expression`. This is a real, currently-unresolved incompatibility between `StanHeaders`/`RcppParallel`'s bundled TBB threading headers and recent Apple Clang — not something fixable by swapping package versions (we tried downgrading `tbb`, which just traded one compile error for another, and forces an unrelated R-version downgrade across the whole conda env — not worth it). `05-run_model.R` sidesteps this on the local branch by using `cmdstanr` instead of `rstan` for `MODEL_TYPE="STAN"` (see [Software prerequisites](#software-prerequisites)); the HPC branch is untouched and keeps using `rstan`, since this Clang-version issue doesn't apply there. If you hit this error, you're likely running the script directly with an R that has `rstan` but not `cmdstanr` set up, or with `HPC` not set to `FALSE` — double check both.

**A script silently uses HPC-only paths and fails with "file not found."** You almost certainly forgot to set `HPC=FALSE` — see [required edits, item 1](#1-the-hpc-environment-variable-defaults-the-wrong-way).

---
