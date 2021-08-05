#' International trade agreements data set
#'
#' A panel data set containing bilateral trade flows between 210 exporters and 262 importers, as
#' well as 50 provision dummies selected at random from a broader data set.
#'
#' @format A data frame with 194092 rows and 59 variables.
#' \itemize{
#'   \item exp: Exporter country (ISO code)
#'   \item imp: Importer country (ISO code)
#'   \item time: Year
#'   \item export: Bilateral trade flows in USD (TODO: check this)
#'   \item id: Agreement ID code
#'   \item agreement: Agreement name
#'   \item fta_wto: TODO: check this
#'   \item EIA_bbf: TODO: check this
#'   \item fta_eia: TODO: check this
#'   \item Other: Provision dummies. TODO: add reference table or link
#' }
#' @source TODO: check this.
"trade"

#' Country ISO Codes
#'
#' An auxiliary data set with basic geographic information about ISO codes included in the \code{trade}
#' data set.
#'
#' @format A data frame with 249 rows and 4 variables.
#' \itemize{
#'   \item iso: Country ISO code.
#'   \item name: Country name.
#'   \item region: Continent.
#'   \item subregion: sub-continental region.

#' }
#' @source https://github.com/lukes/ISO-3166-Countries-with-Regional-Codes.
"countries"
