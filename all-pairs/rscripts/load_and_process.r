# Generates diagrams of the RTD of each pair of ranks.

require(rhdf5)
require(dplyr)

source("analyze_all_pairs.r")

# analyze_data.r $FILENAME
args = commandArgs(trailingOnly=TRUE)
filepath = args[1]
print(args)

# Get all datasets in this file ending with MEDIAN
filestructure = h5ls(filepath)
med_sets = filter(filestructure, grepl('median', name))

res     = apply(med_sets, 1, function(name){print(name["name"])})

path_list   = unlist(strsplit(filepath, "/"))
filename    = path_list[length(path_list)]
name_list   = unlist(strsplit(filename, ".", fixed = TRUE)) 
fwithoutext = head(name_list, -1)
 
print(paste(sep=" ", "Analyze File", filename))
print(fwithoutext)

if(!dir.exists(fwithoutext)) {
  dir.create(fwithoutext)
}

sapply(res, function(name){
  data = h5read(filepath, name);
  kernel_list = unlist(strsplit(name, "_"))
  kernel      = paste(head(kernel_list, -1), collapse=" ")
  print(kernel)
  analyze_all_pairs(data, kernel, fwithoutext);
})


