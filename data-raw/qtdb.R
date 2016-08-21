# Select ECG files from PhysioBank 'qtdb'
allfiles <- list.files("data-raw/qtdb/")
ids <- unlist(lapply(allfiles, function(x) strsplit(x, "-")[[1]][1]))
for(id in ids){

  # Save ECG signal and annotations as object
  file.phys <- paste0("data-raw/qtdb/", id, "-phys.txt")
  ekg <- read.delim(file.phys, sep = "")
  file.q1c <- paste0("data-raw/qtdb/", id, "-q1c.txt")
  atr <- read.delim(file.q1c, sep = "")
  assign(id, list(ekg, atr))
}

# Bundle as data in package
do.call(devtools::use_data, ids)
