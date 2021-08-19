#' International trade agreements data set
#'
#' A panel data set containing bilateral trade flows between 210 exporters and 262 importers between
#' 1964 and 2016. The data set also contains information about trade agreements in force between
#' country pairs, as well as 16 dummies for specific provisions in those agreements (a small selection
#' from a broader data set).
#'
#' @format A data frame with 194,092 rows and 22 variables:
#' \describe{
#'   \item{exp}{Exporter country (ISO 3166 code)}
#'   \item{imp}{Importer country (ISO 3166 code).}
#'   \item{time}{Year.}
#'   \item{export}{Merchandise trade exports in USD.}
#'   \item{id}{Agreement ID code.}
#'   \item{agreement}{Agreement name.}
#'   \item{ad_prov_14}{Anti-dumping actions allowed and with specific provisions for material injury.}
#'   \item{cp_prov_23}{Does the agreement contain provisions that promote transparency?}
#'   \item{tbt_prov_07}{Technical Regulations - Is the use of international standards promoted?}
#'   \item{tbt_prov_33}{Does the agreement go beyond the TBT (Technical Barriers to Trade) Agreement?}
#'   \item{tf_prov_41}{Harmonization and common legal framework}
#'   \item{tf_prov_45}{Issuance of proof of origin}
#'   \item{ser_prov_47}{Does the agreement contain a standstill provision?}
#'   \item{inv_prov_22}{Does the agreement grant Fair and Equitable Treatment (FET)?}
#'   \item{et_prov_38}{Prohibits export-related performance requirements, subject to exemptions.}
#'   \item{ipr_prov_44}{Stipulates that GIs can be registered and protected through a TM system}
#'   \item{env_prov_18}{Does the agreement require states to control ozone-depleting substances?}
#'   \item{ipr_prov_15}{Incorporates/reaffirms all multilateral agreements to which both parties are a
#'       party (general obligation)}
#'   \item{moc_prov_21}{Does the transfer provision explicitly exclude “good faith and non-discriminatory
#'       application of its laws” related to bankruptcy, insolvency or creditor rights protection?}
#'   \item{ste_prov_30}{Does the agreement regulate subsidization to state enterprises?}
#'   \item{lm_prov_10}{Does the agreement include reference to internationally recognized labor standards?}
#'   \item{cp_prov_26}{Does the agreement regulate consumer protection?}
#' }
#' @source Data on international trade flows was obtained from Comtrade.
#' Provision data comes from:
#' Mattoo, A., N. Rocha, M. Ruta (2020). Handbook of deep trade agreements. Washington, DC: World Bank.
"trade"

#' Country ISO Codes
#'
#' An auxiliary data set with basic geographic information about country ISO 3166 codes included in the
#' \code{trade} data set.
#'
#' @format A data frame with 249 rows and 4 variables.
#' \itemize{
#'   \item iso: Country ISO 3166 code.
#'   \item name: Country name.
#'   \item region: Continent.
#'   \item subregion: sub-continental region.

#' }
#' @source The source of the data set is Luke Duncalfe's ISO-3166-Countries-with-Regional-Codes repository
#' on GitHub (\url{https://github.com/lukes/ISO-3166-Countries-with-Regional-Codes#readme}).
"countries"
