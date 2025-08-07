library(R.utils)

# For each sample, drosophila melanogaster (dm.sam) reads were counted using
reads[i] <- countLines("dm.sam")

# The minimum number of reads was identified
minimum <- min(reads)

# Ratios were computed
ratios <- reads / minimum

# For data with SE and PE samples, SE numbers were doubled to compare with PE numbers
