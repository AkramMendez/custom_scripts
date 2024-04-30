library(tidyverse)
library(data.table)
library(caret)


# Loading reference dataset (COX-2 Activity Data) -------------------------
# COX-2 Activity Data
# From Sutherland, O’Brien, and Weaver (2003): A set of 467 cyclooxygenase-2 (COX-2) inhibitors has been assembled from the published work of a single research group, with in vitro activities against human recombinant enzyme expressed as IC50 values ranging from 1 nM to >100 uM (53 compounds have indeterminate IC50 values).
# 
# A set of 255 descriptors (MOE2D and QikProp) were generated. To classify the data, we used a cutoff of 2^{2.5} to determine activity.
# 
# Using data(cox2) exposes three R objects: cox2Descr is a data frame with the descriptor data, cox2IC50 is a numeric vector of IC50 assay values and cox2Class is a factor vector with the activity results.

data(cox2)
