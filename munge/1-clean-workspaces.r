# The $FIL entries in the XML workspaces are incorrect. They use the *.txt file
# extension rather than the *.fcs extension. This prevents flowWorkspace from
# being able to construct a flowSet properly because the wrong filename is
# being used.
FCS_path <- "data/FCS-files"
workspace_path <- "data"
cache_path <- "cache"

# Each of the folders of FCS files corresponds to a patient ID.
# We extract the patient IDs from the folder names. 
# Due to a bug in `list.dirs`, the full paths are returned even if the
# `full.names` argument is set to FALSE.
patient_IDs <- list.dirs(path = FCS_path, full.names = FALSE, recursive = FALSE)
patient_IDs <- sapply(strsplit(patient_IDs, split = "/"), tail, n = 1)

# We copy the XML workspaces to the `cache_path` before cleaning them with the
# `sed` command to fix the references to the FCS files. 
for (patient_ID in patient_IDs) {
  workspace_file <- paste0(patient_ID, ".xml")

  file.copy(from = file.path(workspace_path, workspace_file),
            to = file.path(cache_path, workspace_file))
  
  workspace_file <- file.path(cache_path, workspace_file)
  
  # For the current patient, determine all the FCS files.
  FCS_files <- dir(path = file.path(FCS_path, patient_ID), pattern = ".fcs")
  
  # For each FCS file, call `sed` to fix the XML file.
  # The -i.backup argument backs up the XML file before updating it inline.
  # This is useful in case an error occurs.
  # Example usage:
  # sed -i.bak s/D1NSpre_cells_found.txt/D1NSpre_cells_found.fcs/g 110226.xml
  for (FCS_file in FCS_files) {
    txt_file <- gsub("fcs", "txt", FCS_file, fixed = TRUE)
    sed_command <- paste("s", txt_file, FCS_file, "g", sep = "/")
    sed_command <- paste("sed -i.backup", sed_command, workspace_file)
    system(sed_command)
  }
}

# If all goes well, the XML workspace files have been edited successfully.
# The `sed` command has created *.xml.backup files in `cache_path`.
# We remove those here.
backup_files <- dir(path = cache_path, pattern = "backup", full.names = TRUE)
file.remove(backup_files)
