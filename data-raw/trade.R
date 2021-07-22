## This is the code to prepare 'trade' data set for export. This file will be ignored when the package
# is built from source:

rm(list = ls())
require(haven)

# Trade data (Total trade, every 4 years)
WB_TRADE_DATA <- read_dta("/Users/diego/Documents/OneDrive - London School of Economics/CEP/Summer RA 2020/Files as of June 2021 (for Diego)/Data/temp_trade_only.dta")

# Large provisions data
WB_LARGE   <- read_dta("/Users/diego/Documents/OneDrive - London School of Economics/CEP/Summer RA 2020/Files as of June 2021 (for Diego)/Data/temp_provisions_largedataset_essential_Jan302021.dta")

# Merge:

trade <- merge(WB_TRADE_DATA,WB_LARGE,c("iso1","iso2","year"),all.x=TRUE,all.y=FALSE)
rm(WB_LARGE)
rm(WB_TRADE_DATA)

# Save as .Rdata file, with compression:

usethis::use_data(trade, overwrite = TRUE, compress = "xz")
