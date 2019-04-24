BTSPAS_input <- function(recovery, commercial, 
                         rel.stratum, rec.stratum, 
                         stratum.index, catch.var="CdnCommCt"){

  if(length(recovery[,rel.stratum])==0 | length(recovery[,rec.stratum])==0){
    stop("Input strata variables don't exists in recovery data frame")
  }
  if(length(commercial[,rec.stratum])==0){
    stop("Recovery strata variables don't exists in commercial data frame")
  }  
  #browser()
  # only want recoveries from commercial fishery. If the RecoveryType is anything else
  # then convert the rec.stratum to missing
  recovery[,rec.stratum][!recovery$RecoveryType %in% "Commercial"] <- NA
  
  # convert release stratum to character form and check for missing values
  # if the stratum has an "NA" in the stratum value, it is also considered to be missing
  recovery[, rel.stratum] <- as.character(recovery[, rel.stratum])
  # Check that no missing values in release statum
  if(sum(grepl("NA",recovery[,rel.stratum]) | is.na(recovery[,rel.stratum]))){
    stop("Missing values not allowed in release stratum")
  }
  # Check for no missing values in commercial recovery strata
  if(sum(grepl("NA",commercial[,rec.stratum]) | is.na(commercial[,rec.stratum]))){
    stop("Missing values not allowed in commercial recoveries")
  }
  
  #browser()
  # convert release stratum to character form and set to missing any stratum value that contains NA 
  recovery[, rec.stratum] <- as.character(recovery[, rec.stratum])
  recovery[, rec.stratum][ grepl("NA", recovery[, rec.stratum])] <- NA
  
  # only keep those release strata that appears in statum.index$labels
  recovery <- recovery[ recovery[,rel.stratum] %in% stratum.index$stratum.label,]
  # only keep those recovery strata that appears in statum.index$labels
  recovery <- recovery  [ recovery  [,rec.stratum] %in% stratum.index$stratum.label | is.na(recovery[,rec.stratum]),]
  commercial<-commercial[ commercial[,rec.stratum] %in% stratum.index$stratum.label,]
  
  # add the stratum index values to original data frame for releases and recoveries
  recovery <- merge(recovery, plyr::rename(stratum.index, c("stratum.index"="rel.index")), 
                    by.x=rel.stratum, by.y="stratum.label", all.x=TRUE)
  recovery <- merge(recovery, plyr::rename(stratum.index, c("stratum.index"="rec.index")), 
                    by.x=rec.stratum, by.y="stratum.label", all.x=TRUE)
  
  # compute statistics
  
  # total number of releases by the stratum index
  n1.df <- plyr::ddply(recovery, "rel.index", plyr::summarize, n1=length(rel.index))
  # impute any zeroes for no releases in a stratum. We only impute up to larges stratum index
  max.rel.index <- max(n1.df$rel.index)
  n1.df <- merge(n1.df,  plyr::rename(stratum.index[,"stratum.index", drop=FALSE], c("stratum.index"="rel.index")), all.y=TRUE)
  n1.df <- n1.df[ n1.df$rel.index <= max.rel.index,]
  n1.df$n1[ is.na(n1.df$n1)] <- 0
  n1.df <- n1.df[ order(n1.df$rel.index),]
  #browser()
  # get the cross classification of releases x recoveries. We first need to insert 0's in missing combinations
  m2.large <- plyr::ddply(recovery, "rel.index", function(x){
    m2.df <- plyr::ddply(x, "rec.index", plyr::summarize, m2=length(rec.index))
    # impute any zeroes for no recoveries in a stratum. 
    m2.df <- merge(m2.df,  plyr::rename(stratum.index[,"stratum.index", drop=FALSE], c("stratum.index"="rec.index")),
                   all.y=TRUE)
    m2.df$m2[ is.na(m2.df$m2)] <- 0
    m2.df
  })
  # insert 0 for any release strata with no tags put out
  m2.all.comb <- expand.grid(rel.index= n1.df$rel.index, rec.index=stratum.index$stratum.index)
  m2.large <- merge(m2.large, m2.all.comb, all.y=TRUE)
  m2.large$m2[ is.na(m2.large$m2)] <- 0
  
  # now we cross classify
  m2.full <- as.matrix(xtabs(m2~rel.index+rec.index, data=m2.large))
  
  # rotate each row to the left because you can't recover before releases
  # we put the amount of (shift+1) in the first column
  m2.red <- plyr::aaply(cbind(n1.df$rel.index,m2.full), 1, function(x){
    #browser()
    rot.vec <- c((x[1]+1):length(x), 2:(x[1]+1))
    x <- x[rot.vec]
    x[-length(x)]
  })
  #browser()
  # remove columns at the left that are all zero
  all.zero <- apply(m2.red==0, 2, all)
  remove.right <- rev(cumprod(rev(all.zero)))
  m2.red <- m2.red[, !remove.right]
  
  # remove columns at the left that are all zero
  all.zero <- apply(m2.red==0, 2, all)
  remove.left <- cumprod(all.zero)
  m2.red <- m2.red[, !remove.left]
  
  # Find the commercial recoveries by rec.stratum
  # Merge the stratum index
  # browser()
  commercial <- merge(commercial, plyr::rename(stratum.index, c("stratum.index"="rec.index")), 
                      by.x=rec.stratum, by.y="stratum.label", all.x=TRUE)
  #browser()
  u2.df <- plyr::ddply(commercial, "rec.index", function(x){
    u2=sum(x[,catch.var])
    data.frame(u2=u2)
  })
  # impute zeros for missing recovery strata
  u2.df <-merge(u2.df,  plyr::rename(stratum.index[,"stratum.index", drop=FALSE], c("stratum.index"="rec.index")), all.y=TRUE)
  u2.df$u2[ is.na(u2.df$u2)] <- 0
  Qu2.df <- u2.df[ order(u2.df$rec.index),]
  #browser()
  commercial.catch.df <- u2.df 
  # subtract the recoveries of marked fish
  u2.df.before.correction <- u2.df$u2 - apply(m2.full,2,sum)
  u2.df$u2 <- pmax(0,u2.df.before.correction)
  
  #browser()
  # round everything to integer values
  n1.df <- round(n1.df)
  m2.full <- round(m2.full)
  m2.red  <- round(m2.red)
  commercial.catch.df <- round(commercial.catch.df)
  u2.df.before.correction <- round(u2.df.before.correction)
  u2.df <- round(u2.df)
  
  list(Year=commercial$Year[1],
       recovery=recovery,
       commercial=commercial, 
       stratum.index=stratum.index,
       n1.df=n1.df,
       m2.full=m2.full,
       m2.red=m2.red,
       commercial.catch.df=commercial.catch.df,
       u2.df.before.correction=u2.df.before.correction,
       u2.df=u2.df)
}


# Make the call to BTSPAS
fit.BTSPAS <- function( input.data, prefix="BTSPAS", debug=FALSE, increase.iterations.factor=1){
  # many.more.iterations = increase interations by factor of x
  # takes the input data created earlier and fits the BTSPAS model
  #browser()
  cat("Starting to fit year ", input.data$Year, "\n")
  dir.name <- paste(prefix,"-", input.data$Year, sep="")
  cat("Directory name is ", dir.name, "\n")
  if(file.access(dir.name)!=0){
    dir.create(dir.name, showWarnings=TRUE)
  }
  setwd(dir.name)
  taku.prefix <- paste(prefix,"-",input.data$Year, sep="")
  taku.title  <- paste(taku.prefix," TSPND NP")
  
  
  #write out a csv file with the releases and recoveries in one giant data frame
  #browser()
  n1.df.copy <- input.data$n1.df
  m2.full.copy <- as.data.frame.matrix(input.data$m2.full)
  colnames(m2.full.copy) <- input.data$stratum.index$stratum.label[1:ncol(m2.full.copy)]
  
  out.df <- cbind(n1.df.copy, m2.full.copy)
  u2.copy <- as.data.frame.matrix(t(input.data$u2.df$u2))
  colnames(u2.copy)<- input.data$stratum.index$stratum.label[1:ncol(m2.full.copy)]
  u2.copy$n1 <- NA
  u2.copy$rel.index <- NA
  
  total.cap <- as.data.frame(t(unlist(apply(input.data$m2.full,2,sum) + input.data$u2.df$u2)))
  colnames(total.cap) <- input.data$stratum.index$stratum.label[1:ncol(m2.full.copy)]
  total.cap$n1 <- NA
  total.cap$rel.index <- NA
  
  out.df <- rbind(out.df, u2.copy, total.cap)
  row.names(out.df) <- c(input.data$stratum.index$stratum.label[1:nrow(n1.df.copy)],"untagged","total.captures")
  
  write.csv(out.df, paste(taku.prefix,"-SummaryStatisticsInMatrixForm.csv",sep=""), row.names=TRUE)
  
  
  # what is the strata identification number (statistical week from start of year)?
  taku.sweek <- input.data$stratum.index$stratum.index
  
  # releases, recoveries, and commercial catch
  taku.n1 <- input.data$n1.df$n1
  taku.m2 <- input.data$m2.red  # we want the reduced matrix
  taku.u2 <- input.data$u2.df$u2
  
  # are there any jumps in the abundance?
  taku.jump.after <- NULL    # list sample times after which jump in number occurs
  
  # are there any bad values that need to be adjusted?
  taku.bad.n1     <- c()     # list sample times of bad n1 values
  taku.bad.m2     <- c()     # list sample times of bad m2 values
  taku.bad.u2     <- c()     # list sample times of bad u2 values
  
  # are there any days where the capture probability is fixed in advance?, i.e. because no commercial fishery
  #  taku.logitP.fixed        <- taku.sweek[which(taku.u2 ==0)]
  taku.logitP.fixed        <- taku.sweek[which(total.cap ==0)]
  taku.logitP.fixed.values <- rep(-10, length(taku.logitP.fixed))
  #browser()
  ## Run TSPNDE model
  taku.fit.tspndenp <- TimeStratPetersenNonDiagErrorNP_fit(
    title=      taku.title,
    prefix=     taku.prefix,
    time=       taku.sweek,
    n1=         taku.n1,
    m2=         taku.m2,
    u2=         taku.u2,
    jump.after= taku.jump.after,
    bad.n1=     taku.bad.n1,
    bad.m2=     taku.bad.m2,
    bad.u2=     taku.bad.u2,
    logitP.fixed=taku.logitP.fixed,
    logitP.fixed.values=taku.logitP.fixed.values,
    n.iter=ifelse(debug,1000,20000*increase.iterations.factor), 
    n.burnin=ifelse(debug,100,1000), n.sims=ifelse(debug,30,500),
    debug=FALSE
  )
  # Rename files that were created.
  
  file.rename("data.txt",       paste(taku.prefix,".data.txt",sep=""))
  file.rename("CODAindex.txt",  paste(taku.prefix,".CODAindex.txt",sep=""))
  file.rename("CODAchain1.txt", paste(taku.prefix,".CODAchain1.txt",sep=""))
  file.rename("CODAchain2.txt", paste(taku.prefix,".CODAchain2.txt",sep=""))
  file.rename("CODAchain3.txt", paste(taku.prefix,".CODAchain3.txt",sep=""))
  file.rename("inits1.txt",     paste(taku.prefix,".inits1.txt",sep=""))
  file.rename("inits2.txt",     paste(taku.prefix,".inits2.txt",sep=""))
  file.rename("inits3.txt",     paste(taku.prefix,".inits3.txt",sep=""))
  
  
  #Save the information for later retreival if needed
  save(list=c("taku.fit.tspndenp"), file="taku-fit-tspndenp-saved.Rdata")  # save the results from this run
  setwd("..")# return back
}


# Make the call to BTSPAS allowing for dropout/fallback
fit.BTSPAS.dropout <- function( input.data, prefix="BTSPAS", debug=FALSE, n, dropout){
  # takes the input data created earlier and fits the BTSPAS model
  #browser()
  cat("Starting to fit year ", input.data$Year, "\n")
  dir.name <- paste(prefix,"-", input.data$Year, sep="")
  cat("Directory name is ", dir.name, "\n")
  if(file.access(dir.name)!=0){
    dir.create(dir.name, showWarnings=TRUE)
  }
  setwd(dir.name)
  taku.prefix <- paste(prefix,"-",input.data$Year, sep="")
  taku.title  <- paste(taku.prefix," TSPND NP")
  
  
  #write out a csv file with the releases and recoveries in one giant data frame
  #browser()
  n1.df.copy <- input.data$n1.df
  m2.full.copy <- as.data.frame.matrix(input.data$m2.full)
  colnames(m2.full.copy) <- input.data$stratum.index$stratum.label[1:ncol(m2.full.copy)]
  
  out.df <- cbind(n1.df.copy, m2.full.copy)
  u2.copy <- as.data.frame.matrix(t(input.data$u2.df$u2))
  colnames(u2.copy)<- input.data$stratum.index$stratum.label[1:ncol(m2.full.copy)]
  u2.copy$n1 <- NA
  u2.copy$rel.index <- NA
  
  total.cap <- as.data.frame(t(unlist(apply(input.data$m2.full,2,sum) + input.data$u2.df$u2)))
  colnames(total.cap) <- input.data$stratum.index$stratum.label[1:ncol(m2.full.copy)]
  total.cap$n1 <- NA
  total.cap$rel.index <- NA
  
  out.df <- rbind(out.df, u2.copy, total.cap)
  row.names(out.df) <- c(input.data$stratum.index$stratum.label[1:nrow(n1.df.copy)],"untagged","total.captures")
  
  write.csv(out.df, paste(taku.prefix,"-SummaryStatisticsInMatrixForm.csv",sep=""), row.names=TRUE)
  
  
  # what is the strata identification number (statistical week from start of year)?
  taku.sweek <- input.data$stratum.index$stratum.index
  
  # releases, recoveries, and commercial catch
  taku.n1 <- input.data$n1.df$n1
  taku.m2 <- input.data$m2.red  # we want the reduced matrix
  taku.u2 <- input.data$u2.df$u2
  
  # are there any jumps in the abundance?
  taku.jump.after <- NULL    # list sample times after which jump in number occurs
  
  # are there any bad values that need to be adjusted?
  taku.bad.n1     <- c()     # list sample times of bad n1 values
  taku.bad.m2     <- c()     # list sample times of bad m2 values
  taku.bad.u2     <- c()     # list sample times of bad u2 values
  
  # are there any days where the capture probability is fixed in advance?, i.e. because no commercial fishery
  #taku.logitP.fixed        <- taku.sweek[which(taku.u2 ==0)]
  taku.logitP.fixed        <- taku.sweek[which(total.cap ==0)]
  taku.logitP.fixed.values <- rep(-10, length(taku.logitP.fixed))
  #browser()
  ## Run TSPNDE model
  taku.fit.tspndenp <- TimeStratPetersenNonDiagErrorNPMarkAvail_fit(
    title=      taku.title,
    prefix=     taku.prefix,
    time=       taku.sweek,
    n1=         taku.n1,
    m2=         taku.m2,
    u2=         taku.u2,
    jump.after= taku.jump.after,
    bad.n1=     taku.bad.n1,
    bad.m2=     taku.bad.m2,
    bad.u2=     taku.bad.u2,
    logitP.fixed=taku.logitP.fixed,
    logitP.fixed.values=taku.logitP.fixed.values,
    n.iter=ifelse(debug,1000,30000), n.burnin=ifelse(debug,100,1000), n.sims=ifelse(debug,30,500),
    marked_available_n=n, marked_available_x=n-dropout,
    debug=FALSE
  )
  # Rename files that were created.
  
  file.rename("data.txt",       paste(taku.prefix,".data.txt",sep=""))
  file.rename("CODAindex.txt",  paste(taku.prefix,".CODAindex.txt",sep=""))
  file.rename("CODAchain1.txt", paste(taku.prefix,".CODAchain1.txt",sep=""))
  file.rename("CODAchain2.txt", paste(taku.prefix,".CODAchain2.txt",sep=""))
  file.rename("CODAchain3.txt", paste(taku.prefix,".CODAchain3.txt",sep=""))
  file.rename("inits1.txt",     paste(taku.prefix,".inits1.txt",sep=""))
  file.rename("inits2.txt",     paste(taku.prefix,".inits2.txt",sep=""))
  file.rename("inits3.txt",     paste(taku.prefix,".inits3.txt",sep=""))
  
  
  # Save the information for later retreival if needed
  save(list=c("taku.fit.tspndenp"), file="taku-fit-tspndenp-saved.Rdata")  # save the results from this run
  setwd("..")# return back
}
