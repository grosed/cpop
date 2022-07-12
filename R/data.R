#' Wavenumber Power Spectra data.
#'
#' Data of power spectra of velocity as a function of wavenumber obtained from climate models of the Atlantic Ocean for two different months and
#' two climate scenarios: a present-day one (run 2000) and a future scenario (run 2100). See \insertRef{richards-whitt}{cpop}  
#' 
#' @name wavenumber_spectra
#' 
#' @docType data
#'
#' @keywords datasets
#'
#' @usage data(wavenumber_spectra)
#'
#' @rdname wavenumber-spectra-data
#'
#' @format A dataframe with 5 columns and 247 rows. The first column is wavenumber (cycles per metre). 
#' The other columns are the power spectra values for different months (Feb run 2000, Aug run 2000, Feb run 2100, Aug run 2100).
#' The original data is documented in \insertRef{richards-whitt-data}{cpop}
#' 
#' @examples
#' library(cpop)
#' library(pacman)
#' p_load(tidyr,ggplot2,dplyr)
#'
#' data(wavenumber_spectra)
#'
#' # take logs of variables
#' data <-  wavenumber_spectra %>%  mutate_all(log) %>% rename_all( ~ paste0("log_", .x))
#' head(data)
#'
#' # reproduce figure 4 in "The Impact of Climate Change on Ocean Submesoscale 
#' # Activity" - Richards and Whitt (2021)
#' data %>%
#' gather(variable,log_power_spectra,-log_wavenumber) %>%
#' ggplot(aes(x=log_wavenumber, y=log_power_spectra, colour=variable)) +
#' geom_line() + theme_bw()
#'
#' @references \insertRef{richards-whitt}{cpop}
#' @references \insertRef{richards-whitt-data}{cpop}
#'
NULL












