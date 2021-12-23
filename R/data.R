#' Wavenumber Power Spectra data.
#'
#' Data of power spectra of velocity as a function of wavenumber obtained from climate models of the Atlantic Ocean at two depths and for two different months \insertCite{richards-whitt}{cpop}
#' \insertCite{richards_kelvin_j_2020_4615129}{cpop}. 
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
#' The other columns are the power spectra values for different months (Feb depth 100m, Aug depth 100m, Feb depth 50m, and  Aug depth 50m).
#' The original data can be downloaded from  https://zenodo.org/record/4615129/files/spectra_36hr_05_1.25km_all.mat?download=1 
#' 
#' @examples
#  # reproduce figure 4 in "The Impact of Climate Change on Ocean Submesoscale Activity" - Richards and Whitt (2021) 
#' library(cpop)
#' library(pacman)
#' p_load(tidyr,ggplot2,dplyr)
#' data(wavenumber_spectra)
#' # take logs of variables
#' data <-  wavenumber_spectra %>%  mutate_all(log) %>% rename_all( ~ paste0("log_", .x))
#' head(data)
#' # reproduce figure 4 in Richards and Whitt 
#' data %>%
#' gather(variable,log_power_spectra,-log_wavenumber) %>%
#' ggplot(aes(x=log_wavenumber, y=log_power_spectra, colour=variable)) +
#' geom_line() + 
#' theme_bw()
#' @references \insertRef{richards-whitt}{cpop} \insertRef{richards_kelvin_j_2020_4615129}{cpop}
#'
NULL








