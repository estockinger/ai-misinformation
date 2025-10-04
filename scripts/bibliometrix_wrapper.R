library(ggplot2) 
library(stringr)


topKWs <- function(M, Tag="DE", synonyms=c(), remove.terms=c(), top=10, sep = ";", treatGroups=function(x, ...) {x}) {
  synonyms_indiv <- strsplit(synonyms, ";")
  remove.terms<-remove.terms |> map_if(function(i) any(grepl(i, synonyms_indiv)), function(i) synonyms_indiv[grepl(i, synonyms_indiv)]) %>% unlist %>% unique
  
  topExtKW <- KeywordGrowthMod(M, Tag=Tag, remove.terms=remove.terms, synonyms=synonyms, top=top, cdf=FALSE) %>%
    gather(Tag, Occurrences, -Year) %>%
    group_by(Year) %>% 
    mutate(Extrap = if_else(Year == 2025, Occurrences, 0)) %>%
    filter(Tag != "NA") %>%
    mutate(Frequency=sum(Occurrences)) %>% 
    mutate(Frequency_Prog=sum(Extrap)) %>% 
    group_by(Year,Tag) %>% 
    mutate(Frequency=Occurrences/Frequency) %>%
    mutate(Frequency_Prog=(Occurrences+Extrap)/Frequency_Prog) %>%
    group_by(Tag) %>%
    mutate(CDF = cumsum(Occurrences)) %>% 
    mutate(CDF_Prog = cumsum(Occurrences+Extrap)) %>% 
    ungroup() %>% 
    replace(is.na(.), 0)
}

# quick and dirty work-around to use word boundaries in the native bibliometrix function, which capitalized escaped characters
KeywordGrowthMod <- function(M, Tag = "ID", sep = ";", top = 10, cdf = TRUE, remove.terms = NULL, synonyms = NULL, wordlim=2) {
  i <- which(names(M) == Tag)
  PY <- as.numeric(M$PY)
  Tab <- (strsplit(as.character(M[, i]), sep))
  Y <- rep(PY, lengths(Tab))
  A <- data.frame(Tab = unlist(Tab), Y = Y)
  A$Tab <- trim.leading(A$Tab)
  A <- A[A$Tab != "", ]
  A <- A[!is.na(A$Y), ]
  
  ### remove terms
  terms <- data.frame(Tab = toupper(remove.terms))
  A <- anti_join(A, terms)
  # end of block
  
  ### Merge synonyms in the vector synonyms
  if (length(synonyms) > 0 & is.character(synonyms)) {
    s <- strsplit(toupper(synonyms), ";") |> 
      purrr::map(function (x) { modify_if(.x=x, .p = function(x) {nchar(x) <= wordlim}, .f = function(x) {paste0("\\b",x,"\\b")}) })
    snew <- trimws(unlist(lapply(s, function(l) l[1])))
    sold <- (lapply(s, function(l) trimws(l[-1])))
    for (i in 1:length(s)) {
      A <- A %>%
        mutate(
          Tab = stringi::stri_replace_all_regex(Tab, stringi::stri_replace_all_regex(stringi::stri_replace_all_regex(paste(sold[[i]], collapse = "|", sep = ""), "\\(", "\\\\("), "\\)", "\\\\)"), snew[i])
        )
    }
  }
  # end of block
  Ymin <- min(A$Y)
  Ymax <- max(A$Y)
  Year <- Ymin:Ymax
  Tab <- names(sort(table(A$Tab), decreasing = TRUE))[1:top]
  
  words <- matrix(0, length(Year), top + 1)
  words <- data.frame(words)
  names(words) <- c("Year", Tab)
  words[, 1] <- Year
  for (j in 1:length(Tab)) {
    word <- (table(A[A$Tab %in% Tab[j], 2]))
    words[, j + 1] <- trim.years(word, Year, cdf)
  }
  return(words)
}

trim.years <- function(w, Year, cdf) {
  Y <- as.numeric(names(w))
  W <- matrix(0, length(Year), 1)
  
  for (i in 1:length(Year)) {
    if (Y[1] == Year[i] & length(Y) > 0) {
      W[i, 1] <- w[1]
      Y <- Y[-1]
      w <- w[-1]
    }
  }
  if (isTRUE(cdf)) W <- cumsum(W)
  names(W) <- Year
  W <- data.frame(W)
  return(W)
}