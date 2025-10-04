library(countrycode)
library(dplyr)

uniCountryAffil <- function(data) {
  uniCountryDf <- data[, c("AU1_CO", "AU1_UN", "AU_CO", "AU_UN")] |> 
    mutate(across(everything(), ~gsub("UNIV ", "",gsub(" UNIVERSITY", "", gsub("UNIVERSITY OF ", "", str_extract(.x, "[^;]+$"))))))
  uniCountryDf <- data.frame(Map(c,uniCountryDf[, 1:2], uniCountryDf[, 3:4])) |> unique()
  colnames(uniCountryDf) <- c("country", "uni")
  
  uniCountryDf[!( # manually correct wrong affiliations
    grepl("CORRESPONDING", uniCountryDf$uni) |
      (uniCountryDf$country != "DENMARK" & uniCountryDf$uni == "SOUTHERN DENMARK") |
      (uniCountryDf$country != "UNITED KINGDOM" & uniCountryDf$uni %in% c("EXETER", "EAST LONDON", "OXFORD")) |
      (uniCountryDf$country != "SWITZERLAND" & uniCountryDf$uni == "ZURICH") |
      (uniCountryDf$country != "CHINA" & uniCountryDf$uni == "NOTTINGHAM NINGBO CHINA") |
      (uniCountryDf$country != "AUSTRALIA" & uniCountryDf$uni == "MELBOURNE")
  ), ] |>
    na.omit()
}

countryCountryCollabs <- function(data){
  uniCountryDf <- uniCountryAffil(data)
  collabDf <- data[, c("AU1_UN", "AU_UN")] |> 
    mutate(across(everything(), ~gsub("UNIV ", "",gsub(" UNIVERSITY", "", gsub("UNIVERSITY OF ", "", str_extract(.x, "[^;]+$")))))) |> 
    na.omit()
  
  collabDf <- collabDf[!(grepl("CORRESPONDING", collabDf$AU1_UN) | grepl("CORRESPONDING", collabDf$AU_UN)), ] |>
    left_join(uniCountryDf, by = c("AU1_UN" = "uni")) |> 
    left_join(uniCountryDf, by = c("AU_UN" = "uni")) |>
    group_by(country.x, country.y) |>
    summarize(value = n(), .groups = "drop") |>
    na.omit()
  colnames(collabDf) <- c("from", "to", "value")
  collabDf
}

countryCollabStats <- function(countryCountryDf, replace_smallest_n=25) {
  countryCollabs <- full_join(
    countryCountryDf |> summarise(value = sum(value), .by=from), 
    countryCountryDf |> summarise(value = sum(value), .by=to), 
    by=(c("from"="to")))
  colnames(countryCollabs) <- c("country", "from", "to")
  
  countryCollabs$overall <- countryCollabs$from + countryCollabs$to
  countryCollabs[is.na(countryCollabs)] <- 0
  countryCollabs$continent <- countrycode(sourcevar = countryCollabs$country, origin = "country.name", destination = "continent")
  
  # replace smallest countries with continent values
  smallestCountries <- countryCollabs |> arrange(overall) |> slice_head(n=replace_smallest_n) |> pull(country)
  countryCollabs <- countryCollabs |> mutate(chordEntry = if_else(country %in% smallestCountries, continent, country))
  countryCollabs
}