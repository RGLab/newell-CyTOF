library(flowWorkspace)

FCS_path <- "data/FCS-files"
workspace_path <- "data"
cache_path <- "cache"

# Each of the folders of FCS files corresponds to a patient ID.
# We extract the patient IDs from the folder names. 
# Due to a bug in `list.dirs`, the full paths are returned even if the
# `full.names` argument is set to FALSE.
patient_IDs <- list.dirs(path = FCS_path, full.names = FALSE, recursive = FALSE)
patient_IDs <- sapply(strsplit(patient_IDs, split = "/"), tail, n = 1)

workspace_files <- dir(path = cache_path, pattern = "xml", full.names = TRUE)

message("Parsing workspaces...")
# Parses workspaces and constructs a flowSet
patient_ID <- patient_IDs[1]
fs_list <- sapply(patient_IDs, function(patient_ID) {
  message("Patient: ", patient_ID)
  
  workspace_file <- paste0(file.path(cache_path, patient_ID), ".xml")
  fcs_path <- file.path(FCS_path, patient_ID)
  
  ws <- openWorkspace(workspace_file)
  gating_set <- parseWorkspace(ws, name = 1, path = fcs_path, isNcdf = FALSE)
  closeWorkspace(ws)
  
  getData(gating_set)
}, simplify = FALSE)