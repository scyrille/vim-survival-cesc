
#'@description Data preprocessing

library(here)

source(here::here("R","utils.R"))

# Import data -------------------------------------------------------------

raids <- read_processed("raids")
tcga <- read_processed("tcga")

# Combine data ------------------------------------------------------------

combined <- combine_data(raids, tcga)
saveRDS(combined, here::here("data","processed","combined",
                             "combined_processed_data.RDS"))

