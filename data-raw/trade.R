## This is the code to prepare 'trade' data set for export. This file will be ignored when the package
# is built from source:

rm(list = ls())
require(haven)
require(stringi)

# Trade data (Total trade, every 4 years)
WB_TRADE_DATA <- read_dta("C:/Users/LocalAdmin/Dropbox/Machine Learning and International Trade/Data/temp_trade_only.dta")

# Large provisions data
WB_LARGE <- read_dta("C:/Users/LocalAdmin/Dropbox/Machine Learning and International Trade/Data/temp_provisions_largedataset_essential_Jan302021.dta")

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

# Finally, to make the data set lighter, we drop most of the provision data and keep just the
# six selected provisions plus 10 others chosen at random:

provisions <- names(trade[, -(1:9)])
selected <- c("ad_prov_14", "cp_prov_23", "tbt_prov_07", "tbt_prov_33", "tf_prov_41", "tf_prov_45")
randomly_selected <- c("env_prov_18", "et_prov_38", "inv_prov_22",
              "ipr_prov_15", "ipr_prov_44", "lm_prov_10", "moc_prov_21",
              "ser_prov_47", "ste_prov_30", "cp_prov_26")
# Previously randomly_selected was selected with this, but since
# no seed was set, I choose this here.
#provisions <- provisions[!(provisions %in% selected)]
#selected <- c(selected, sample(provisions, size = 10))#
selected <- c(selected, randomly_selected)
trade <- trade[, c("exp", "imp", "time", "export", "id", "agreement", selected)]
trade$agreement <- iconv(trade$agreement, "utf-8", "ASCII", sub="")

# Save as .Rdata file, with compression to stay below CRAN's 5Mb size limit:
usethis::use_data(trade, compress = "xz", overwrite = TRUE)

# We may also want to include a complementary data set with information about specific countries. This one
# provides ISO code, name, region and subregion:

countries <- read.csv("~/Documents/OneDrive - London School of Economics/CEP/Summer RA 2020/ISO-3166-Countries-with-Regional-Codes-9.0/all/all.csv")
countries <- countries[, c(3, 1, 6, 7)]
names(countries)[1] <- c("iso")

# We ensure all country names use only ASCII characters (otherwise R CMD check will throw a warning):

countries$name <- stringi::stri_escape_unicode(countries$name)

# We could also add an OECD variable (TODO).

# Finally, we save it:

usethis::use_data(countries, compress = "xz", overwrite = TRUE)

