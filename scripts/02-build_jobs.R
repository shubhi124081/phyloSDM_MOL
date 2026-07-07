# Purpose of script: Build bash script that will run a specified script (specified)
# in config

# Set up ----
HPC <- Sys.getenv("HPC")
if (HPC == "FALSE") {
    root <- "~/clim_risk_phylosdm"
} else {
    root <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm"
}

hpc_root_filepath <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm"
# Filepaths & scripts
scripts_directory <- file.path(root, "scripts")
data_directory <- file.path(root, "data")
job_directory <- file.path(root, "jobs")
log_directory <- file.path(root, "log")

##################### Read Config  #####################
# Configuration parameters for data generation
config <- yaml::read_yaml(file.path(scripts_directory, "01-config.yaml"),
    as.named.list = TRUE
)

data <- config$data
script <- config$script
tmb_specs <- config$tmb_specs
model_specs <- config$model_specs
hpc_specs <- config$hpc_specs
buildLocal <- data$buildLocal

# Extract model parameters - error if missing
if (is.null(model_specs$model_type)) {
    stop("ERROR: model_specs$model_type is required in config but is NULL")
}
if (is.null(model_specs$model_name) || length(model_specs$model_name) == 0) {
    stop("ERROR: model_specs$model_name is required in config but is NULL or empty")
}
model_type <- model_specs$model_type
model_name <- model_specs$model_name[[1]]

# Extract job parameters from config
cluster_name_config <- data$cluster_name
repno_config <- data$repno
focal_sp_config <- data$focal_sp
nrep <- data$nrep
raw_data_name <- data$raw_data
exp_id <- data$exp_id
exp_root <- data$exp_root
shortname <- exp_id

##################### Load Cluster Data #####################
# Load the cluster mapping data
spList_path <- file.path(
    root, "raw_data", raw_data_name, "spList.Rdata"
)
contents <- load(spList_path)
# This loads: spList

##################### Expand Clusters #####################
# 1. Cluster name expansion
# spList is a named list where names are cluster names and values are species vectors
all_cluster_names <- names(spList)

if (is.null(cluster_name_config) ||
    toupper(as.character(cluster_name_config)) == "ALL") {
    # Get all unique cluster names
    clusters_to_run <- all_cluster_names
    message(sprintf("Expanding to all %d clusters", length(clusters_to_run)))
} else {
    # Use the specified cluster(s)
    clusters_to_run <- cluster_name_config
    # Validate that specified clusters exist
    invalid_clusters <- setdiff(clusters_to_run, all_cluster_names)
    if (length(invalid_clusters) > 0) {
        stop(sprintf(
            "Invalid cluster name(s): %s. Available clusters: %s",
            paste(invalid_clusters, collapse = ", "),
            paste(head(all_cluster_names, 10), collapse = ", ")
        ))
    }
}

##################### Expand Repno #####################
# 2. Repno expansion
if (toupper(as.character(repno_config)) == "ALL") {
    repnos_to_run <- seq_len(nrep)
    message(sprintf("Expanding to all %d reps", nrep))
} else {
    repno_val <- as.integer(repno_config)
    # Validate repno is not greater than nrep

    if (repno_val > nrep) {
        stop(sprintf(
            "repno (%d) cannot be greater than nrep (%d)",
            repno_val, nrep
        ))
    }
    if (repno_val < 1) {
        stop(sprintf("repno (%d) must be at least 1", repno_val))
    }
    repnos_to_run <- repno_val
    message(sprintf("Number of reps to run %d", length(repnos_to_run)))
}

##################### Expand Focal Species #####################
# 3. Focal species expansion
# Get all species within the specified clusters from spList
valid_species <- unlist(spList[clusters_to_run], use.names = FALSE)

focal_sp_upper <- toupper(as.character(focal_sp_config))

if (focal_sp_upper == "ALL") {
    # "ALL" means cluster-level job - pass "ALL" to script, one job per cluster
    # The downstream script will handle all species within the cluster
    focal_sps_to_run <- "ALL"
    message(sprintf(
        "Cluster-level mode: creating one job per cluster (script handles species)"
    ))
} else if (focal_sp_upper == "RUN_ALL_SEP") {
    # "run_all_sep" means separate job for each species
    focal_sps_to_run <- unique(valid_species)
    message(sprintf(
        "Expanding to all %d focal species across %d cluster(s)",
        length(focal_sps_to_run), length(clusters_to_run)
    ))
} else {
    # Use the specified focal species
    focal_sps_to_run <- focal_sp_config
    # Validate that focal species exist in the specified clusters
    invalid_species <- setdiff(focal_sps_to_run, valid_species)
    if (length(invalid_species) > 0) {
        stop(sprintf(
            "Invalid focal species: %s. Species not found in specified cluster(s).",
            paste(invalid_species, collapse = ", ")
        ))
    }
}

##################### Build Experiment Files Grid #####################
# Create a grid of all combinations to run
# Each job will have: cluster_name, repno, focal_sp

if (length(focal_sps_to_run) == 1 && focal_sps_to_run == "ALL") {
    # Cluster-level mode: one job per cluster with focal_sp = "ALL"
    focal_cluster_map <- data.frame(
        focal_sp = "ALL",
        cluster_name = clusters_to_run
    )
} else {
    # Species-level mode: build mapping from spList
    focal_cluster_map <- do.call(rbind, lapply(clusters_to_run, function(cl) {
        species_in_cluster <- spList[[cl]]
        species_to_use <- intersect(species_in_cluster, focal_sps_to_run)
        if (length(species_to_use) > 0) {
            data.frame(focal_sp = species_to_use, cluster_name = cl)
        } else {
            NULL
        }
    }))
}
focal_cluster_map <- unique(focal_cluster_map)

# Expand grid with repnos
expfiles <- merge(
    focal_cluster_map,
    data.frame(repno = repnos_to_run),
    by = NULL # Cartesian product
)

# Add all metadata columns
expfiles$dataset <- raw_data_name
expfiles$exp_root <- exp_root
expfiles$exp_id <- exp_id
expfiles$nrep <- nrep

# Reorder columns for clarity
expfiles <- expfiles[, c("exp_root", "exp_id", "dataset", "cluster_name", "focal_sp", "repno", "nrep")]

message(sprintf(
    "Built %d job combinations: %d cluster(s) x %d species x %d rep(s)",
    nrow(expfiles),
    length(clusters_to_run),
    length(focal_sps_to_run),
    length(repnos_to_run)
))

Rscript_prefix <- "Rscript "
# Generally "Rscript " will work, but in case you need to point to another
# install of Rscript, change this

header <- "
module load miniconda
conda activate phylo-sdms_wf_clone"
shebang <- "#!/bin/bash"
dsq <- hpc_specs$dsq

# Two ways to submit - job array or job on loop
# way one
if (dsq) {
    modules <- "module load miniconda; conda activate phylo-sdms_wf_clone; "
    Rscript_prefix <- "Rscript "

    # Use the script specified in config
    script_name <- script$script_name
    command <- paste0(
        Rscript_prefix,
        file.path(hpc_root_filepath, "scripts", script_name)
    )
    # Output file
    script1 <- gsub("\\.R", "", script_name)
    exp_name1 <- shortname
    filename <- paste0("job_array_", exp_name1, "_for_", script1, ".txt")
    output_file <- file.path(job_directory, filename)
    nloop <- nrow(expfiles)

    # Open a connection to the file
    file_conn <- file(output_file, open = "wt")

    for (o in seq_len(nloop)) {
        # Pass all job parameters: exp_root, exp_id, dataset, cluster_name, focal_sp, repno, nrep, model_type, model_name
        command_fin <- paste(
            modules, command,
            expfiles[o, "exp_root"],
            expfiles[o, "exp_id"],
            expfiles[o, "dataset"],
            expfiles[o, "cluster_name"],
            expfiles[o, "focal_sp"],
            expfiles[o, "repno"],
            expfiles[o, "nrep"],
            model_type,
            model_name, "\n"
        )
        cat(sprintf(command_fin), file = file_conn)
    }

    # Close the file connection
    close(file_conn)

    print(sprintf(
        "To submit, load dsq and create a job file with name %s",
        filename
    ))

    system(sprintf("job_name=%s", filename))
} else {
    # Use the script specified in config
    script_name <- script$script_name

    if (buildLocal) {
        # If built local but run on the HPC, then write the hpc filepath into the
        # job script
        command <- paste0(
            Rscript_prefix,
            file.path(hpc_root_filepath, "scripts", script_name),
            " "
        )
    } else {
        command <- paste0(
            Rscript_prefix,
            file.path(scripts_directory, script_name),
            " "
        )
    }

    prefix <- "#SBATCH --"

    for (o in seq_len(nrow(expfiles))) {
        # Create descriptive job name using cluster and species
        jobname_main <- paste0(
            shortname, "_",
            expfiles[o, "cluster_name"], "_",
            expfiles[o, "focal_sp"], "_rep",
            expfiles[o, "repno"], ".sh"
        )
        job_file <- file(file.path(job_directory, jobname_main), "w")

        mid <- paste0(
            prefix, "job-name=", jobname_main, "\n",
            prefix, "partition=", hpc_specs$PARTITION, "\n",
            prefix, "out=", "log/slurm-%A_%a.out", "\n",
            prefix, "time=", hpc_specs$TIME, "\n",
            prefix, "nodes=", hpc_specs$NODES, "\n",
            prefix, "ntasks=", hpc_specs$NTASKS, "\n",
            prefix, "cpus-per-task=", hpc_specs$CPUSPERTASK, "\n",
            prefix, "mem=", hpc_specs$MEM, "\n",
            prefix, "mail-user=", hpc_specs$MAIL_USER, "\n",
            prefix, "mail-type=", hpc_specs$MAIL_TYPE
        )

        # Pass all job parameters: exp_root, exp_id, dataset, cluster_name, focal_sp, repno, nrep, model_type, model_name
        command_fin <- paste(
            command,
            expfiles[o, "exp_root"],
            expfiles[o, "exp_id"],
            expfiles[o, "dataset"],
            expfiles[o, "cluster_name"],
            expfiles[o, "focal_sp"],
            expfiles[o, "repno"],
            expfiles[o, "nrep"],
            model_type,
            model_name
        )

        writeLines(
            c(shebang, mid, "\n", header, command_fin),
            file.path(job_directory, jobname_main)
        )

        if (hpc_specs$AUTO_SUBMIT) {
            system(paste("sbatch", file.path(job_directory, jobname_main)))
        }

        if (!hpc_specs$AUTO_SUBMIT) {
            print(paste0(
                "To submit, run sbatch ", file.path(job_directory, jobname_main)
            ))
        }
    }
}

# Save the job manifest for reference
write.csv(expfiles,
    file = file.path(log_directory, paste0(shortname, "_jobs.csv")),
    row.names = TRUE
)
message(sprintf(
    "Job manifest saved to: %s",
    file.path(log_directory, paste0(shortname, "_jobs.csv"))
))
