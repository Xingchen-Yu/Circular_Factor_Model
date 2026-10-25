######The only difference between readKH2 in this script and readKH from "wnominate" pacakge lies in line 61 to 64 
####in which the two non-working functions commented out are replaced by trim.leading and trim.trailing defined in line 5 and 8.

# returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x)

# returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)

readKH2<-function (file, dtl = NULL, yea = c(1, 2, 3), nay = c(4, 5, 
                                                               6), missing = c(7, 8, 9), notInLegis = 0, desc = NULL, debug = FALSE) 
{
  cat("Attempting to read file in Keith Poole/Howard Rosenthal (KH) format.\n")
  warnLevel <- options()$warn
  options(warn = -1)
  data <- try(readLines(con = file), silent = TRUE)
  if (inherits(data, "try-error")) {
    cat(paste("Could not read", file, "\n"))
    return(invisible(NULL))
  }
  options(warn = warnLevel)
  cat("Attempting to create roll call object\n")
  voteData <- substring(data, 37)
  n <- length(voteData)
  m <- nchar(voteData)[1]
  rollCallMatrix <- matrix(NA, n, m)
  for (i in 1:n) {
    asdf<-gsub(" ","",voteData[i])
    rollCallMatrix[i, ] <- as.numeric(unlist(strsplit(asdf, 
                                                      split = character(0))))
  }
  rm(voteData)
  if (!is.null(desc)) 
    cat(paste(desc, "\n"))
  cat(paste(n, "legislators and", m, "roll calls\n"))
  cat("Frequency counts for vote types:\n")
  tab <- table(rollCallMatrix, exclude = NULL)
  print(tab)
  icpsrLegis <- as.numeric(substring(data, 4, 8))
  party <- as.numeric(substring(data, 21, 23))
  partyfunc <- function(x) {
    party <- partycodes$party[match(x, partycodes$code)]
    party[party == "Democrat"] <- "D"
    party[party == "Republican"] <- "R"
    party[party == "Independent"] <- "Indep"
    party
  }
  partyName <- partyfunc(party)
  statename <- function(x) {
    state.info$state[match(x, state.info$icpsr)]
  }
  state <- as.numeric(substring(data, 9, 10))
  KHstateName <- substring(data, 13, 20)
  stateName <- statename(state)
  stateAbb <- datasets::state.abb[match(stateName, datasets::state.name)]
  stateAbb[grep(KHstateName, pattern = "^USA")] <- "USA"
  cd <- as.numeric(substring(data, 11, 12))
  cdChar <- as.character(cd)
  cdChar[cd == 0] <- ""
  lnames <- substring(data, 26, 36)
  for (i in 1:n) {
    #lnames[i] <- strip.trailing.space(lnames[i])
    lnames[i] <- trim.leading(lnames[i])
    lnames[i] <- trim.trailing(lnames[i])
    #lnames[i] <- strip.after.comma(lnames[i])
  }
  legisId <- paste(lnames, " (", partyName, " ", stateAbb, 
                   "-", cdChar, ")", sep = "")
  legisId <- gsub(x = legisId, pattern = "-)", replacement = ")")
  if (any(duplicated(legisId))) {
    dups <- duplicated(legisId)
    legisId[dups] <- paste(legisId[dups], icpsrLegis[dups])
  }
  legis.data <- data.frame(state = stateAbb, icpsrState = state, 
                           cd = cd, icpsrLegis = icpsrLegis, party = partyName, 
                           partyCode = party)
  dimnames(legis.data)[[1]] <- legisId
  vote.data <- NULL
  if (!is.null(dtl)) {
    vote.data <- dtlParser(dtl, debug = debug)
  }
  rc <- rollcall(data = rollCallMatrix, yea = yea, nay = nay, 
                 missing = missing, notInLegis = notInLegis, legis.names = legisId, 
                 legis.data = legis.data, vote.data = vote.data, desc = desc, 
                 source = file)
  rc
}