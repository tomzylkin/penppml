## This is the code to prepare 'trade' data set for export. This file will be ignored when the package
# is built from source:

rm(list = ls())
require(haven)

# Trade data (Total trade, every 4 years)
WB_TRADE_DATA <- read_dta("/Users/diego/Documents/OneDrive - London School of Economics/CEP/Summer RA 2020/Files as of June 2021 (for Diego)/Data/temp_trade_only.dta")

# Large provisions data
WB_LARGE <- read_dta("/Users/diego/Documents/OneDrive - London School of Economics/CEP/Summer RA 2020/Files as of June 2021 (for Diego)/Data/temp_provisions_largedataset_essential_Jan302021.dta")


# Merge:

trade <- merge(WB_TRADE_DATA,WB_LARGE, c("iso1", "iso2", "year"), all.x = TRUE, all.y = FALSE)
rm(WB_LARGE)
rm(WB_TRADE_DATA)

# Light cleaning:

trade[is.na(trade)] <- 0     # Remove NAs.
trade <- trade[, 1:314]      # Last few variables are not provision dummies.

# We drop observations for agreements no longer in effect (we do not have provisions data for these):

excluded <- ((rowSums(data.matrix(trade[, -1:-9]))) == 0) & (trade$fta_eia == 1)
trade <- trade[!excluded, ]

# Transform ID variables into factors (and rename them so they are more informative):

trade$iso1  <- factor(trade$iso1)
trade$iso2  <- factor(trade$iso2)
trade$year <- ordered(trade$year)
colnames(trade)[1:3] <- c("exp", "imp", "time")

# Finally, to facilitate testing, we drop most of the provision data and keep just the iceberg provisions:

iceberg <- read.csv("/Users/diego/Documents/OneDrive - London School of Economics/CEP/Summer RA 2020/Files as of June 2021 (for Diego)/Data/icberg_provisions.csv", header = FALSE)
trade <- trade[, c("exp", "imp", "time", "export", "id", "agreement", iceberg$V1[iceberg$V2 == 1])]

# Save as .Rdata file, with compression to stay below CRAN's 5Mb size limit:

usethis::use_data(trade, compress = "xz", overwrite = TRUE)
