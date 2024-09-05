setwd("/Users/thibaultbret/dbcan-output")

library(dplyr)
source("/Users/thibaultbret/Documents/Work/dbcan4 output family cleaner.R")

### Metagenome profiles ####
#===========================#

## Reading in the general data
filelist = Sys.glob(paste(getwd(), "/*dbcan-sub.hmm.out", sep = ""))

# Reading in all data files
subs = lapply(filelist, read.table, header = TRUE, sep = "\t", stringsAsFactors=FALSE, row.names = NULL)

substrate_profiles = NULL
for (i in 1:length(filelist)){
  current_profile = get_substrates(subs[[i]])
  if (nrow(current_profile) == 0) {
    cat("No data found for file:", filelist[i], "\n")
    next
  }
  current_profile$bin = sub('.*output/(.*)\\.dbcan.*', '\\1', filelist[i])
  substrate_profiles = rbind(substrate_profiles, current_profile)
}

write.csv(substrate_profiles, "/Users/thibaultbret/substrate_profiles.csv")
saveRDS(substrate_profiles, "/Users/thibaultbret/substrate_profiles.RDS")