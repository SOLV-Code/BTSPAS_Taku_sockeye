# Inseason Estimates for Taku River
#
#  After reading in the data and doing various merges,
#  select the stat weeks for which you want the BTSPAS to provide estimates
#  on a FW and HW basis.
# 
#  It will create a series of directories in the current workspace that will
#  accumulate as you run the code each week
#  This code will compute the Petersen, the Full Week (FW) and Half week (HW) stratified
#  Petersen using BTSPAS
#



# Data files will be provided on a weekly basis from the tagging crews and DFO
# commercial catch
#
# release data
#   - when fish are released with the following variable names
#        Year, TagID,  ReleaseDate, ReleaseStatWeek (starts on Sunday)

# recapture data from DFO
#   - which tags are recovered in COMMERCIAL catch only with the following names
#        Year, TagID, RecoveryDate, RecoveryStatWeek (starts on Sunday), RecoveryType
#     only those records with RecoveryType="Commerical" will be used

# commercial catch from DFO'
#   - commercial catch EXCLUDING recoveries of tagged fish with the following names
#        Year, RecoveryDate, RecoveryStatWeek CatchWithTags, RecoveryType
#     only those records with RecoveryType="Commerical" will be used
#     the commercial catch should INCLUDE the count of tagged fish recovered
#     **** It is assumed that the recovery date matches a commerical opening. For example,
#     if a tag is returned after the opening is closed, it is assume to have occurred during
#     the opening (which is usally in the first half of the week) ****
#     THIS IS IMPORTANT TO GET THE HALF WEEK ANALYSIS TO WORK. See the checks later in the code

##### BE SURE TO download the latest version of BTSPAS from
##### the GitHub site at https://github.com/cschwarz-stat-sfu-ca/BTSPAS
##### using
#####     devtools::install_github("cschwarz-stat-sfu-ca/BTSPAS", 
#####                              dependencies = TRUE,
#####                              build_vignettes = TRUE)
##### This could take up to 20 minutes, so be patient. (The vignettes take a long time
##### to compile.

library(BTSPAS)
library(ggplot2)
library(lubridate)

source("../myfunctions.r")

#---------------------
# read in the data and make sure that variable names match
# the list above
release <- read.csv(file.path("InSeasonData","release_data.csv"),
                    header=TRUE, as.is=TRUE, strip.white=TRUE)
release$ReleaseDate <- lubridate::ymd(release$ReleaseDate)
dim(release)
release <- release[ !is.na(release$ReleaseDate),]
dim(release)
head(release)

recap <- read.csv(file.path("InSeasonData","recovery_data.csv"),
                  header=TRUE, as.is=TRUE, strip.white=TRUE)
recap$RecoveryDate <- lubridate::ymd(recap$RecoveryDate)
recap$RecoveryType <- "Commercial"
head(recap)

catch <- read.csv(file.path("InSeasonData","catch_data.csv"),
                  header=TRUE, as.is=TRUE, strip.white=TRUE)
catch$Date <- lubridate::ymd(catch$Date)
catch <- plyr::rename(catch, c("Date"="RecoveryDate",
                               "StatWeek"="RecoveryStatWeek",
                               "CdnCommCt"="CatchWithTags"))
catch$RecoveryType="Commercial"
head(catch)


# some basic statistics
xtabs(~ReleaseStatWeek,  data=release,exclude=NULL, na.action=na.pass)
xtabs(~RecoveryStatWeek, data=recap,  exclude=NULL, na.action=na.pass)
xtabs(CatchWithTags~RecoveryStatWeek, data=catch,  exclude=NULL, na.action=na.pass)

#------------------------------------------
# Figure out the half week strata
# We will divide Sunday -> Wednesday as the first half= Commercial Opening
#                Thursday -> Saturday as the second half
#
# First get all the dates in the study
strata <- rbind( 
  plyr::rename(release, c("ReleaseStatWeek"="StatWeek",
                          "ReleaseDate"="Date"))[,c("Year","Date","StatWeek")],
  plyr::rename(recap,   c("RecoveryStatWeek"="StatWeek",
                          "RecoveryDate"="Date"))[,c("Year","Date","StatWeek")],
  plyr::rename(catch,   c("RecoveryStatWeek"="StatWeek",
                          "RecoveryDate"="Date"))[,c("Year","Date","StatWeek")])
strata <- unique(strata[,c("Year","Date","StatWeek")])

strata <- strata[ order(strata$Date),]
# Get the day of the week,
strata$dow <- lubridate::wday(strata$Date,week_start=7)

# Days 1-4 (Sunday -> Wed) = first half of week = Commercial Opening
# Days 5-7 (Thurs  -> Sat) = last half of week
strata$HalfStatWeek <- strata$StatWeek + .1 +
  .1*strata$dow %in% c(5,6,7)
head(strata)

# mergback the strata halfweek back to release, recovery,
# and commercial catch

dim(release)
release <- merge(release, 
                 plyr::rename(strata, 
                              c("Date"="ReleaseDate",
                                "HalfStatWeek"="ReleaseHalfStatWeek"))
                 [,c("ReleaseDate","ReleaseHalfStatWeek")],
                 all.x=TRUE)
dim(release)

dim(recap)
recap <- merge(recap, 
               plyr::rename(strata, 
                            c("Date"="RecoveryDate",
                              "HalfStatWeek"="RecoveryHalfStatWeek"))
               [,c("RecoveryDate","RecoveryHalfStatWeek")],
               all.x=TRUE)
dim(recap)

# but a recap that occurs in the second half of the week is assumed to have occured
# in the commerical opening. You need to change the data as needed for the half week
# analysis and is likely to be year specific
xtabs(~RecoveryHalfStatWeek, data=recap)

select <- (recap$RecoveryHalfStatWeek %% 1) > .15  # recoveries in second half of week
sum(select)
xtabs(~RecoveryHalfStatWeek, data=recap[select,])
recap$RecoveryHalfStatWeek[select] <- recap$RecoveryHalfStatWeek[select] - .1 # shif back to the opening

xtabs(~RecoveryHalfStatWeek, data=recap)


dim(catch)
catch <- merge(catch, 
               plyr::rename(strata, 
                            c("Date"="RecoveryDate",
                              "HalfStatWeek"="RecoveryHalfStatWeek"))
               [,c("RecoveryDate","RecoveryHalfStatWeek")],
               all.x=TRUE)
dim(catch)

#----------------------------------------

# merge the recaptures with the releases
# check that recapture tag numbers match
setdiff(recap$TagID, release$TagID)

dim(release)
relrecap <- merge(release, recap, all.x=TRUE)
dim(relrecap)

#----------------------------------------
#----------------------------------------
#----------------------------------------
#----------------------------------------

# Generate the inseason estimate
# Give the list the StatWeeks that should be included in the estimate

xtabs(~ReleaseStatWeek,  data=relrecap)
xtabs(~RecoveryStatWeek, data=relrecap)
xtabs(CatchWithTags~RecoveryStatWeek, data=catch) 


#-------------
# Full Week BTSPAS analysis
# Define the stratum variable as 1 = first stat week, 2=second stat week etc
# THIS IS WHERE YOU SELECT THE STAT WEEKS FROM YOUR DATA SET TO MAKE THE ESTIMATES
# REFER TO THE XTABS() JUST ABOVE FOR SOME HELP IN MAKING THE DECISION

fw.stat.weeks <- 23:28   # stat weeks with releases and recoveries to  be included

fw.stratum.index <- data.frame(stratum.index=1:length(fw.stat.weeks),
                               stratum.label=as.character(fw.stat.weeks),
                               stringsAsFactors=FALSE)
fw.stratum.index


# get the data necessary for the call to BTSPAS
fw.data <- BTSPAS_input(relrecap, catch, "ReleaseStatWeek", "RecoveryStatWeek",
                        fw.stratum.index, catch.var="CatchWithTags")

# fit the BTSPAS model
fw.prefix <- paste("Taku-FW-Inseason-W",round(min(fw.stat.weeks)),
                   "-W",round(max(fw.stat.weeks)),"-",sep="")
fit.BTSPAS(fw.data,prefix=fw.prefix)

# fit the BTSPAS model with fall back (say n=50, x=12)
fw.prefix.dropout <- paste("Taku-FW-Inseason-W",round(min(fw.stat.weeks)),
                           "-W",round(max(fw.stat.weeks)),"-fallback-",sep="")
fit.BTSPAS.dropout(fw.data,prefix=fw.prefix.dropout, n=50, dropout=12)


#-------------
# Half Week BTSPAS analysis
# Define the stratum variable as 1 = first stat week, 2=second stat week etc
hw.stat.weeks <- sort(as.vector(outer(fw.stat.weeks, c(.1,.2), "+")))  # releases and recoveries 
hw.stat.weeks

hw.stratum.index <- data.frame(stratum.index=1:length(hw.stat.weeks),
                               stratum.label=as.character(hw.stat.weeks),
                               stringsAsFactors=FALSE)
hw.stratum.index


# get the data necessary for the call to BTSPAS
hw.data <- BTSPAS_input(relrecap, catch, "ReleaseHalfStatWeek", "RecoveryHalfStatWeek",
                        hw.stratum.index, catch.var="CatchWithTags")

# fit the BTSPAS model
hw.prefix <- gsub("FW","HW",fw.prefix)
fit.BTSPAS(hw.data,prefix=hw.prefix)


# fit the BTSPAS model with fall back (say n=50, x=12)
hw.prefix.dropout <- gsub("FW","HW",fw.prefix.dropout)
fit.BTSPAS.dropout(hw.data,prefix=hw.prefix.dropout, n=50, dropout=12)


#------------------------------------
#------------------------------------
#-----------------------------------
# Make a table of the estimates from the various sets of weeks etc
##### Extract the outputs

## Extract the results from the various fits
file.names <-dir()
# Extract the directories with the fits
file.names.fits<- file.names[grepl(paste("^Taku-"), file.names)]
file.names.fits

# make a pdf file of the fitted curves
pdf(paste("Inseason-all-fits.pdf",sep=""))
plyr::l_ply(file.names.fits, function(x){
  cat("Extracting final plot from ", x, "\n")
  load(file.path(x, "taku-fit-tspndenp-saved.Rdata"))
  tryCatch(plot(taku.fit.tspndenp$plots$fit.plot),
           error=function(cond){
             message("Unable to draw plot - convergence error?- skipped")
             return(NULL)
           })
})
dev.off()

# Extract all of the estimates of the total run size
run.size <- plyr::ldply(file.names.fits, function(x){
  cat("Extracting total run size from ", x, "\n")
  load(file.path(x, "taku-fit-tspndenp-saved.Rdata"))
  Ntot <- taku.fit.tspndenp$summary["Ntot",]
  if(is.null(Ntot))Ntot <- rep(NA,9) # if model didn't converge properly
  #browser()
  Ntot <- as.data.frame(t(Ntot))
  Ntot[,1:7] <- round(Ntot[,1:7])
  Ntot$file=x
  Ntot
})
run.size
write.csv(run.size,
          file=paste("inseason-run.size.csv",sep=""))

## Extract the Petersen estimators
#
# Extract all of the estimates of the total run size
run.pet.size <- plyr::ldply(file.names.fits, function(x){
  cat("Extracting Petersen from ",x,"\n")
  load(file.path(x, "taku-fit-tspndenp-saved.Rdata"))
  Year <- as.numeric(substring(x, 2+regexpr('--',x,fixed=TRUE)))
  #browser()
  Ntot.pp <- taku.fit.tspndenp$PP$using.all.data
  if(is.null(Ntot.pp$N.se))Ntot.pp <- data.frame(N.est=NA, N.se=NA)
  # see if this included fall back
  Ntot.pp.fallback <- taku.fit.tspndenp$PP$using.all.data.fallback
  if(is.null(Ntot.pp.fallback$N.se))Ntot.pp.fallback <- data.frame(N.est=NA, N.se=NA)
  
  c(Ntot.pp.est=round(Ntot.pp$N.est), Ntot.pp.se=round(Ntot.pp$N.se), 
    Ntot.pp.fallback.est=round(Ntot.pp.fallback$N.est), Ntot.pp.fallback.se=round(Ntot.pp.fallback$N.se),
    file=x)
})
run.pet.size
write.csv(run.pet.size,
          file="Inseason-PP.run.size.csv")