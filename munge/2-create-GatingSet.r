library(ProjectTemplate)
load.project()

# The FCS files for each XML workspace have different number of stains.
# To handle this, we open a workspace and then separately parse the
# workspace for each FCS file in the FCS path.
# We then extract the flowFrame for each FCS file and add it to fs_list.
FCS_path <- "data/FCS-files"
workspace_path <- "data"
cache_path <- "cache"
fs_list <- list()

# Patient 110226
fs_110226 <- list()
gslist_110226 <- list()
patient_ID <- "110226"

fcs_path <- file.path(FCS_path, patient_ID)
fcs_files <- dir(path = fcs_path, pattern = "\\.fcs$")

workspace_file <- paste0(file.path(cache_path, patient_ID), ".xml")
ws <- openWorkspace(workspace_file)

for (fcs_file in fcs_files) {
  gslist_110226[[fcs_file]] <- parseWorkspace(ws, name = 1, path = fcs_path,
                                          isNcdf = FALSE, subset = fcs_file)
  fs_110226[[fcs_file]] <- getData(gslist_110226[[fcs_file]][[1]])
}

closeWorkspace(ws)

# The first two FCS files have slightly different channel names than the rest of
# the FCS files.  Specifically, the first two files are appended with the letter
# "i", while the channels for the remaining FCS files are appended with the
# letter "d". I'm unsure why.
# These differing channel names prevent flowCore from combining the flowFrames
# into a single flowSet. We replace the channel names for the first two
# flowFrames with the channel names for the other flowFrames.
fs_110226[[1]] <- copy_flowframe_channels(fs_110226[[1]], fs_110226[[3]])
fs_110226[[2]] <- copy_flowframe_channels(fs_110226[[2]], fs_110226[[3]])

fs_110226 <- flowSet(fs_110226)

# Updates pData and phenoData
pData(fs_110226)$PTID <- patient_ID
varM <- varMetadata(phenoData(fs_110226))
varM[-1,] <- rownames(varM)[-1]
varMetadata(phenoData(fs_110226)) <- varM


# Patient 110317

# There are no problems with channel names for patient 110317. So there is no
# need to update the channel names like we did above.

fs_110317 <- list()
gslist_110317 <- list()
patient_ID <- "110317"

fcs_path <- file.path(FCS_path, patient_ID)
fcs_files <- dir(path = fcs_path, pattern = "\\.fcs$")

workspace_file <- paste0(file.path(cache_path, patient_ID), ".xml")
ws <- openWorkspace(workspace_file)

for (fcs_file in fcs_files) {
  gslist_110317[[fcs_file]] <- parseWorkspace(ws, name = 1, path = fcs_path, isNcdf = FALSE,
                               subset = fcs_file)
  fs_110317[[fcs_file]] <- getData(gslist_110317[[fcs_file]][[1]])
}

closeWorkspace(ws)

# The files with 'a' for appended need to be combined into a single flowFrame.
# To do this, we remove the letter 'a' from each filename and then combine the
# flowFrames having a common filename into a single flowFrame.

# We first sort the filenames so that the ordering of duplicates helps us
# identify which files are the same after we remove the 'a'.
filenames_sorted <- sort(names(fs_110317))
filenames_stripped <- gsub("a", "", filenames_sorted)
which_duplicates <- which(duplicated(filenames_stripped))

# Now that we have identified the duplicates, we combine each pair of duplicates
# into a single flowFrame. The only difference between the two flowFrames other
# than their filenames are their expression matrices. We simply 'rbind' both
# expression matrices and assign them to the original file. Then, we remove the
# "appended" file.
for (dup in which_duplicates) {
  original_file <- filenames_sorted[dup - 1]
  appended_file <- filenames_sorted[dup]

  # First, we rbind the expression matrices.
  fs_110317[[original_file]]@exprs <- rbind(
                                            fs_110317[[original_file]]@exprs,
                                            fs_110317[[appended_file]]@exprs
                                            )
  # Next, we remove the appended file from the flowSet.
  fs_110317[[appended_file]] <- NULL
}

fs_110317 <- flowSet(fs_110317)

# Updates pData and phenoData
pData(fs_110317)$PTID <- patient_ID
varM <- varMetadata(phenoData(fs_110317))
varM[-1,] <- rownames(varM)[-1]
varMetadata(phenoData(fs_110317)) <- varM


# NOTE: I attempted to add enriched samples to 110317 flowSet, but errors
# occurred because the XML workspaces did not have the samples.

# Constructs GatingSets for each patient
gs_110226 <- GatingSet(fs_110226)
gs_110317 <- GatingSet(fs_110317)


# For both patients, we extract the nodes from the list of GatingSets parsed
# from the workspaces. These nodes are common for all GatingSets within the
# list.
# The 'root' note is removed.

# Extracts gates for patient 110226
# For 110226, we update the gates's channel names, similar to above.
node_populations <- setdiff(getNodes(gslist_110226[[1]][[1]]), "root")

for (node in node_populations) {
  node_gates <- list()
  for (fcs_file in sampleNames(fs_110226)) {
    node_gates[[fcs_file]] <- getGate(gslist_110226[[fcs_file]][[1]], node)

    # Updates channel names in gates
    gate_boundaries <- node_gates[[fcs_file]]@boundaries
    colnames(gate_boundaries) <- gsub("Di", "Dd", colnames(gate_boundaries))
    node_gates[[fcs_file]]@boundaries <- gate_boundaries

    gate_parameters <- parameters(node_gates[[fcs_file]])
    gate_parameters <- gsub("Di", "Dd", gate_parameters)
    names(gate_parameters) <- gsub("Di", "Dd", names(gate_parameters))
    parameters(node_gates[[fcs_file]]) <- gate_parameters
  }
  node_filterList <- filterList(node_gates)
  add(gs_110226, node_filterList)
  recompute(gs_110226)
}

# Extracts gates for patient 110317
node_populations <- setdiff(getNodes(gslist_110317[[1]][[1]]), "root")
for (node in node_populations) {
  node_gates <- list()
  for (fcs_file in sampleNames(fs_110317)) {
    node_gates[[fcs_file]] <- getGate(gslist_110317[[fcs_file]][[1]], node)
  }
  node_filterList <- filterList(node_gates)
  add(gs_110317, node_filterList)
  recompute(gs_110317)
}

# Archives GatingSets on blackrhino
save_gs(gs_110226, "/loc/no-backup/ramey/newell-CyTOF/gs_110226")
save_gs(gs_110317, "/loc/no-backup/ramey/newell-CyTOF/gs_110317")

