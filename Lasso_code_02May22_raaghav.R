library(devtools)
devtools::install_github("tomzylkin/penppml") #install packages
library(penppml) #library
library(readxl)
#install.packages("matrixStats")
library(matrixStats)

load("../Raaghav/FInal_LASSO_provision.Rda")

str(FInal_LASSO_provision)
colnames(FInal_LASSO_provision)

FInal_LASSO_provision$iso3_d.x <- as.factor(FInal_LASSO_provision$iso3_d.x) #turn char to factors
class(FInal_LASSO_provision$iso3_d.x)

FInal_LASSO_provision$iso3_o.x <- as.factor(FInal_LASSO_provision$iso3_o.x) #turn char to factors
class(FInal_LASSO_provision$iso3_o.x)


FInal_LASSO_provision$year.x <- factor(FInal_LASSO_provision$year.x, ordered = TRUE,   #turn into original factor
                                       levels = c("2000","2004", "2008", "2012", "2016"))
class(FInal_LASSO_provision$year.x)

FInal_LASSO_provision$value <- as.numeric(FInal_LASSO_provision$value)
FInal_LASSO_provision$value[FInal_LASSO_provision$value < 0] <- 0
class(FInal_LASSO_provision$value)

Missing_Value_Check <- data.frame(colSums(is.na(x=FInal_LASSO_provision))) #missing value check


str(FInal_LASSO_provision)

FInal_LASSO_provision <- data.frame(FInal_LASSO_provision)

FInal_LASSO_provision_stata <- FInal_LASSO_provision
colnames(FInal_LASSO_provision_stata)[3:8] <- c("iso3_d_x",  "partner_country", "iso3_o_x", "fdi",
"direction", "year_x")
write_dta(FInal_LASSO_provision_stata, "../../stata/FinalLassoRaaghav.dta")

lambdas <- c(0.05, 0.025, 0.01, 0.005, 0.0025, 0.001, 0.0005, 0.00025, 0.0001, 0)
lambdas <- c(0.025, 0.001, 0)
lambdas <- c(0.000005, 0.000002, 0.000001, 0.0000005, 0.0000001)
lambdas <- c(0.05, 0.0025, 0.001, 0.0005, 0.00025, 0.0001, 0.000005, 0.000002, 0.000001, 0.0000005, 0.0000001, 0)

reg1 <- mlfitppml(data = FInal_LASSO_provision,
                  dep = "value",
                  indep = 27:331,
                  fixed = list(c("iso3_d.x", "year.x"),
                               c("iso3_o.x", "year.x"),
                               c("iso3_d.x", "iso3_o.x")),
                  penalty = "lasso",
                  lambdas = lambdas, post=FALSE, tol=1e-5)

hdreg <- hdfeppml(data = FInal_LASSO_provision,
         dep = "value", indep="FTA_DUMMY", fixed = list(c("iso3_d.x", "year.x"),
                                                        c("iso3_o.x", "year.x"),
                                                        c("iso3_d.x", "iso3_o.x")))


s_initials <-  unique(unlist(lapply(stringr::str_split(names(reg1$beta[,3]), pattern="_"), "[[", 1)))

res_mat <- data.frame("Provision class" = s_initials, "Number selected"=NA, "Fraction selected"=NA)

i_s <- 1
for(s in s_initials) {
  s_i <- paste("^",s,sep="")
  res_mat[which(res_mat$Provision.class==s), 2] <- length(grep(s_i,names(which(reg1$beta[,3] != 0))))
  res_mat[which(res_mat$Provision.class==s), 3] <- length(grep(s_i,names(which(reg1$beta[,3] != 0))))/length(grep(s_i,colnames(FInal_LASSO_provision[,27:331])))
}

View(res_mat)

reg2 <- mlfitppml(data = FInal_LASSO_provision,
                  dep = "value",
                  indep = 27:331,
                  fixed = list(c("iso3_d.x", "year.x"),
                               c("iso3_o.x", "year.x"),
                               c("iso3_d.x", "iso3_o.x")),
                  penalty = "lasso",
                  method = "plugin",
                  cluster = c("iso3_d.x", "iso3_o.x"))

results <- data.frame(b_pre = reg2$beta_pre, b = reg2$beta, se = 0)
results$se[!is.na(reg9$beta)] <- reg2$ses
results

### Bootstrap LASSO

library(logr)

# Create temp file location
tmp <- file.path("BootstrapLasso.log")

# Open log
lf <- log_open(tmp)

# Send message to log
log_print("Bootstrap Lasso")

# Print data to log
set.seed(125)
bs_sample_size <- 20000
selected_vars <- NULL
for(bs_i in 1:250) {
  tryCatch({
    log_print(paste("MC rep.: ",bs_i))
    bs_sample_ind <- sort(sample(1:dim(FInal_LASSO_provision)[1], bs_sample_size, replace=FALSE))
    data_bs <- FInal_LASSO_provision[bs_sample_ind,]
    reg2 <- mlfitppml(data = data_bs,
                      dep = "value",
                      indep = 27:331,
                      fixed = list(c("iso3_d.x", "year.x"),
                                   c("iso3_o.x", "year.x"),
                                   c("iso3_d.x", "iso3_o.x")),
                      penalty = "lasso",
                      method = "plugin",
                      cluster = c("iso3_d.x", "iso3_o.x"))
    selected_vars <- c(selected_vars, which(reg2$beta != 0))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
log_print(selected_vars)

log_close()

##### Cross validation #####
data_cv <- FInal_LASSO_provision
fes_cv <- list()
inter_cv <- list(c("iso3_d.x", "year.x"),
                 c("iso3_o.x", "year.x"),
                 c("iso3_d.x", "iso3_o.x"))
for (i in seq_along(inter_cv)) {
  fes_cv[[paste(names(data_cv[, inter_cv[[i]]]), collapse = "_")]] <- interaction(data_cv[, inter_cv[[i]]])
}

mat_fes_cv <- matrix(unlist(fes_cv), ncol=length(fes_cv))
data_fes_cv <- data.frame(data_cv, mat_fes_cv)

obs_before <- dim(data_cv)[1]

for(fe_ind in 1:length(fes_cv)){
  fe_name_temp <- paste("X",fe_ind,sep="")
  print(fe_name_temp)
  #print("Nr obs before")
  #print(dim(data_fes)[1])
  data_sum <- data_fes_cv %>% dplyr::group_by(!!rlang::sym(fe_name_temp)) %>% dplyr::summarise(sum_dep = sum(value))
  data_var <- data_fes_cv %>% dplyr::group_by(!!rlang::sym(fe_name_temp)) %>% dplyr::summarise(var_dep = var(value))
  data_temp <- dplyr::left_join(data_fes_cv, data_sum, by=fe_name_temp)
  data_temp <- dplyr::left_join(data_temp, data_var, by=fe_name_temp)
  # Include observations where sum and variance unequal zero and var not NA, i.e. at least two observations in group
  incl_obs <- which(data_temp$sum_dep!=0 & !is.na(data_temp$var_dep) & data_temp$var_dep != 0)
  data_fes_cv <- data_temp[incl_obs,]
  fes_cv <- lapply(fes_cv, "[", incl_obs)
  print(length(fes_cv[[fe_ind]]))
  data_fes_cv <- data_fes_cv %>% dplyr::select(-any_of(c("sum_dep", "var_dep", fe_name_temp)))
}
data <- data_fes_cv
obs_before - dim(data)[1] # This many observations were dropped

data <- data[!is.na(data$id),]
data[data$id=="NO Agrement",20] <- "0"
data$id <- as.numeric(as.character(data$id))

id <- unique(data[, 20])
nfolds <- 10
unique_ids <- data.frame(id = id, fold = sample(1:nfolds, size = length(id), replace = TRUE))

cross_ids <- merge(data[, 20, drop = FALSE],
                   unique_ids, by = "id", all.x = TRUE)

# Tried out different versions of lambda here:
# Uncomment and select the one mentioned in the comment to reproduce
# lambdas <- c(0.005, 0.0005, 0.00005, 0.000005)
# lambdas <- c(0.000001, 0.0000005, 0.0000001)
# lambdas <- c(0.05, 0.000001, 0.0000005, 0.0000001) #minimal
# #lambdas <- c(0.05, 0.0025, 0.001, 0.0005, 0.00025, 0.0001, 0.000005, 0.000002, 0.000001, 0.0000005, 0.0000001, 0)
# lambdas <- c(0.03, 0.01, 0.009, 0.007)
lambdas <- c(seq(0.04,0.01,by=-0.01), 0.008, 0.006)

reg3 <- mlfitppml(data = data,
                  dep = "value",
                  indep = 27:331,
                  fixed = list(c("iso3_d.x", "year.x"),
                               c("iso3_o.x", "year.x"),
                               c("iso3_d.x", "iso3_o.x")),
                  penalty = "lasso",
                  lambdas = lambdas,
                  xval = TRUE,
                  IDs =  cross_ids$fold, post=FALSE)

# lambdas <- 0.05, 0.000001, 0.0000005, 0.0000001
reg3_1 <- reg3
#lambdas <- c(0.005, 0.0005, 0.00005, 0.000005)
reg3_2 <- reg3
# 5e-03 5e-04 5e-05 5e-06
reg3_3 <- reg3

# For this choice of lambdas we select some variables:
#lambdas <- c(seq(0.04,0.01,by=-0.01), 0.008, 0.006)
reg3_4 <- reg3
# ror_prov_03 sps_prov_06  tf_prov_20,  182         220         285

