
# MAKES the qPCR analysis of the raw data (preliminary).
# 
# It gathers all replicae for the different groups from QuantumStudio3 file csv file
# The analysis takes the means and SEM, along with CI for the specified alpha
# 
# The result sheet from QuantumStudio 3 should be exported with tabular separation 
# and quoted text.
# 
# Output files: 
# CTstats - stats on the CT parameter 
# Qstats - stats on the Q parameter (quantity according to the dilution factor)
# 
# @AUTHOR: Maximo Sanchez-Aragon
# @EMAIL: maxsanara@gmail.com
#

args = commandArgs(trailingOnly=TRUE) 

# PARAMS -------------
# experimental folders
wd <- args[1] 
#wd <- "/home/maxkavian" # arg 1

data.fld <- paste(wd,args[2],sep="/") 
#data.fld <- paste(wd,"qPCR_data",sep="/")               # arg 2

results.fld <- paste(wd,args[3],sep="/") 
#results.fld <- paste(wd,"qPCR_analysis",sep="/")        # arg 3

data.exp.fld <- args[4] 
#data.exp.fld <- c("20230705")                           # arg 4

# experimental filename
#data.exp.file <- c("Max_Act88f_rep1_2023-07-05_080147_Results.csv")
#data.exp.file <- c("Max_Act88f_2023-07-05_080147_Results.csv")  # arg 5
data.exp.file <- args[5] #c("Max_RFeSp1_2023-07-09-080800_Results.csv")


# categories ranks that map categs vector and categs.names vector. These come from the plate reader itself
sl <- args[6] #c("a1:a12, b1:d4, b5:d8, e1:g4, e5:g8")  # arg 6
categs.arg <- args[7] #"st.curve,A1_4,A5_8,D1_4,D5_8"   # arg 7
categs.names.arg <- args[8] #c("Standard Curve,tubGS>shRFeSp1 mef(-),tubGS>shRFeSp1 mef(+),tubGS/w[1118] mef(-),tubGS/w[1118] mef(+)")  # arg 8

# OTHER params
stcurve.pos <- as.numeric(args[9]) #as.numeric("1")  # arg 9
skip <- as.numeric(args[10])       #as.numeric("49")    # arg 10, lines in th csv that will be skipped by default, CURRENTLY DEPRECATED *******
alpha <- as.numeric(args[11])      #as.numeric("0.05") # arg 11, alpha for confidence intervals
wgraph <- as.numeric(args[12])     #as.numeric(".75") # arg 12 
wh.arg <- args[13]     #"2,6"             # arg 13
delout <- as.logical(args[14])     # as.logical("TRUE") # arg 14, should outliers be removed?
delhighsd <- as.logical(args[15])  #as.logical("TRUE") # arg 15 should high SD points be removed?



# excluding columns from QuatumStudio3 format
excluding.fields.arg <- args[16] #"1,3,7,8,10,11,13,14,25,30,31,32,33,29"  # arg 16, fields excluded from data... 27 HIGHSD # 28 OUTLIERS
# factoring columns in the file (i.e. make them into factors)
factorings.arg <- args[17] #"Target.Name,Task,Sample.Name,Well.Position"  # arg 17, columns that need refactoring for R compatibility




# param processing
wh <- as.numeric(strsplit(wh.arg,",")[[1]]) # width and height for the graphical output
categs <- strsplit(categs.arg,",")[[1]]
categs.names <- strsplit(categs.names.arg,",")[[1]]
basename<- strsplit(as.character(data.exp.file),".", fixed=T)[[1]][1]
filepath <- paste(data.fld,data.exp.fld,data.exp.file,sep="/")
filepath.res <- paste(results.fld,data.exp.fld,basename,sep="/")
excluding.fields <- as.integer(strsplit(excluding.fields.arg,",")[[1]]) # 27 HIGHSD # 28 OUTLIERS
factorings <- as.character(strsplit(factorings.arg,",")[[1]])

# SKIP lines automatic selection 
skip.selector.string <- "Well Position"
skip <- as.integer(system(paste("cat ",filepath," | grep -n \"",skip.selector.string,"\" | cut -d: -f1",sep=""),intern=TRUE)) -1


# FUNCTIONS ---------
make.rank <- function(s, sep=":"){
  
  # It makes a rank using : as separator in a typical excel notation
  # the result is a character string of the interaction of letters and numbers
  
  s2 <- strsplit(s,sep)[[1]]
  s2ri <- sort(match(toupper(substr(trimws(s2),1,1)), LETTERS))
  s2ci <- sort(as.integer(substr(trimws(s2),2,nchar(trimws(s2)))))

  ris <- LETTERS[c(s2ri[1]:s2ri[2])]
  cis <- c(s2ci[1]:s2ci[2])

  defaultW <- getOption("warn") 
  options(warn = -1) 
  
  f<- levels(interaction(factor(ris),factor(cis), sep=""))
  
  options(warn = defaultW)
  
  f
}
make.rank.list <- function(sl, sepr=",", ...){
  
  # It makes a list of ranks, each corresponding to one category in the original
  # data set
  
  s2 <- trimws(strsplit(sl,sepr)[[1]])
  r <- list()
  for(i in 1:length(s2)){
    r[[i]] <- make.rank(s2[i], ...)
  }
  r
}
make.data <- function(path, categs, categs.ranks, categs.names="", skip=49, excluding.fields=c(), factorings = c(), delete.outliers=TRUE, delete.highsd=FALSE ){
  sl <- categs.ranks
  
  data <- read.table(path,header=TRUE,sep="\t",skip=skip)
  #data <- subset(data,select=-excluding.fields )
  
  for(j in 1:length(factorings))
    data[[factorings[j]]] <- factor(data[[factorings[j]]])

  
  a <- make.rank.list(sl)
  
  categs.cols <- c()
  categs.cols.names <- c()
  for(i in 1:length(categs)){
    categs.cols[match(a[[i]], data[["Well.Position"]] )] <- categs[i]
    categs.cols.names[match(a[[i]], data[["Well.Position"]] )] <- categs.names[i]
  }
  data <- cbind(data, categs.cols,categs.cols.names)
  data$categs.cols <- factor(data$categs.cols)
  data$categs.cols.names <- factor(data$categs.cols.names)
  
  if(delete.highsd){
    if(length(which(names(data) == "HIGHSD") ) > 0)  
    data <- subset(data,HIGHSD!="Y")
  }
  
  if(delete.outliers){
    if(length(which(names(data) == "OUTLIERRG") ) > 0)
    data <- subset(data,OUTLIERRG!="Y")
  }
  
  data
}
SEM.bound <- function(x, q=0.05){
  qt(q, length(x)-1)*SEM.func(x)
}
SEM.func <- function(x){
  sqrt(var(x) / (length(x)-1) )
}
sd <- function(x){
  sqrt(var(x) / (length(x)-1))
}



shapirop <- function(x){
  r <- shapiro.test(x)
  r$p.value
}
make.CTstats <- function(data, alpha=.05, stcurve.pos=1){
  
  # makes the aggregates for the actual measurements, i.e. excluding the standard
  # curve measurements
  datag <- subset(data, categs.cols != categs[stcurve.pos])
  tab.mean <- aggregate(CT ~ categs.cols, datag, mean, na.action = na.omit)
  tab.var <- aggregate(CT ~ categs.cols, datag, var, na.action = na.omit)
  tab.sd <- aggregate(CT ~ categs.cols, datag, sd, na.action = na.omit)
  tab.sem <- aggregate(CT ~ categs.cols, datag, SEM.func, na.action = na.omit)
  tab.semlow <- aggregate(CT ~ categs.cols, datag, SEM.bound, q=alpha/2 ,na.action = na.omit)
  tab.semhigh <- aggregate(CT ~ categs.cols, datag, SEM.bound, q=1-alpha/2 ,na.action = na.omit)
  tab.n <- aggregate(CT ~ categs.cols, datag, length ,na.action = na.omit)
  #tab.shapirop <- aggregate(CT ~ categs.cols, datag, shapirop, na.action=na.action)
  
  names(tab.mean)[2] <- "CT.mean"
  names(tab.var)[2] <- "CT.var"
  names(tab.sd)[2] <- "CT.sd"
  names(tab.sem)[2] <- "CT.SEM"
  names(tab.semlow)[2] <- "CT.SEMlow"
  names(tab.semhigh)[2] <- "CT.SEMhigh"
  names(tab.n)[2] <- "CT.n"
  #names(tab.shapirop)[2] <- "shapiro.p"
  
  tab <- merge(tab.mean,tab.var, by="categs.cols")
  tab <- merge(tab, tab.sd, by="categs.cols")
  tab <- merge(tab, tab.sem, by="categs.cols")
  tab <- merge(tab, tab.semlow, by="categs.cols")
  tab <- merge(tab, tab.semhigh, by="categs.cols")
  tab <- merge(tab, tab.n, by="categs.cols")
  #tab <- merge(tab,tab.shapirop, by="categs.cols")
  
  tab.names <- data.frame(categs.cols = categs, categs.cols.names=categs.names)
  tab.names <-  subset(tab.names, categs.cols != categs[stcurve.pos])
  
  tab <- merge(tab, tab.names, by="categs.cols")
  
  
  attr(tab, "alpha") <- alpha
  
  tab
}
make.Qstats <- function(data, alpha=.05, stcurve.pos=1){
  
  # makes the aggregates for the actual measurements, i.e. excluding the standard
  # curve measurements
  datag <- subset(data, categs.cols != categs[stcurve.pos])
  tab.mean <- aggregate(Quantity ~ categs.cols, datag, mean, na.action = na.omit)
  tab.var <- aggregate(Quantity ~ categs.cols, datag, var, na.action = na.omit)
  tab.sd <- aggregate(Quantity ~ categs.cols, datag, sd, na.action = na.omit)
  tab.sem <- aggregate(Quantity ~ categs.cols, datag, SEM.func, na.action = na.omit)
  tab.semlow <- aggregate(Quantity ~ categs.cols, datag, SEM.bound, q=alpha/2 ,na.action = na.omit)
  tab.semhigh <- aggregate(Quantity ~ categs.cols, datag, SEM.bound, q=1-alpha/2 ,na.action = na.omit)
  tab.n <- aggregate(Quantity ~ categs.cols, datag, length ,na.action = na.omit)
  #tab.shapirop <- aggregate(Quantity ~ categs.cols, datag, shapirop, na.action=na.action)
  
  names(tab.mean)[2] <- "Quantity.mean"
  names(tab.var)[2] <- "Quantity.var"
  names(tab.sd)[2] <- "Quantity.sd"
  names(tab.sem)[2] <- "Quantity.SEM"
  names(tab.semlow)[2] <- "Quantity.SEMlow"
  names(tab.semhigh)[2] <- "Quantity.SEMhigh"
  names(tab.n)[2] <- "Quantity.n"
  #names(tab.shapirop)[2] <- "shapiro.p"
  
  tab <- merge(tab.mean,tab.var, by="categs.cols")
  tab <- merge(tab, tab.sd, by="categs.cols")
  tab <- merge(tab, tab.sem, by="categs.cols")
  tab <- merge(tab, tab.semlow, by="categs.cols")
  tab <- merge(tab, tab.semhigh, by="categs.cols")
  tab <- merge(tab, tab.n, by="categs.cols")
  #tab <- merge(tab,tab.shapirop, by="categs.cols")
  
  tab.names <- data.frame(categs.cols = categs, categs.cols.names=categs.names)
  tab.names <-  subset(tab.names, categs.cols != categs[stcurve.pos])
  
  tab <- merge(tab, tab.names, by="categs.cols")
  
  attr(tab, "alpha") <- alpha
  
  tab
}

# graphical
draw.point <- function(x,y,ydown,yup,w, draw.dot=T, draw.l=T, lty.m="dotted", lwd.m=1, invert=F,rectcoords=NA,  ...){
  
  if(!invert){
    lines(c(x,x),c(ydown,yup), lty=lty.m, lwd = lwd.m, ...)
    if(draw.l){
      lines(c(x-w/2,x+w/2),c(ydown,ydown), lwd=2, ...)
      lines(c(x-w/2,x+w/2),c(yup,yup), lwd=2,  ...)
    }
  }else{
    lines(c(ydown,yup), c(y,y), lty=lty.m, lwd = lwd.m, ...)
    if(draw.l){
      lines(c(ydown,ydown), c(y-w/2,y+w/2), lwd=2, ...)
      lines(c(yup,yup), c(y-w/2,y+w/2), lwd=2,  ...)
    }
  }
  
  if( (!is.na(rectcoords[1])) && (!is.na(rectcoords[2])) && (!invert) ){
    rect(x-w/2, rectcoords[1], x+w/2, rectcoords[2], border=NA, ...)
  }
  
  if(draw.dot){
    points(x, y, pch=21, bg="white", lwd=1, cex=.75, ...)
    #lines(c(x-w/2,x+w/2), c(y,y), col="white", lwd=1, cex=.65)
  }
  
}
make.plot <- function(dt, filepath="", selcolx="categs.cols.names", pad=0.8, wgrap=1, wh=c(2,6), pdf.flag=T, yylab="", make.blank=FALSE, ...){
  selcolx <- "categs.cols.names"
  #dt <- CTstats
  #pad <- 0.8
  #wgraph <- 1
  
  lvs <- dt[[selcolx]]
  x <- 1:length(lvs)
  yylim <- c(min(dt[,c(2)]-dt[,c(7)]), max(dt[,c(2)]+dt[,c(7)]))
  
  if(pdf.flag){
    pdf(file=filepath, wh[1],wh[2])
  }else{
    svg(file=filepath, wh[1],wh[2])
  }
  
  # gaphical
  mar.red <- 1
  # changes the likelihood to switch to scientific notation
  options(scipen = 5)
  #set the distances of the axes labels to the plot
  par(mgp=c(2.4,1,.4))
  # sets the margins
  par(mar=c(12,5-mar.red,2.5-mar.red,2.5-mar.red)+0.1) 
  #par(fig = c(0,1,0,1))
  # sets the size of the character
  par(cex=.9)
  
  plot(NA, xlim=c(min(x)-pad, max(x)+pad), ylim = yylim, type="n", 
       ylab=yylab, xlab="", xaxt="n", ...)
  
  labs <- c()
  for(i in 1:nrow(dt)){
    draw.point(i,dt[i,2],dt[i,2]+dt[i,6],dt[i,2]+dt[i,7],wgraph, draw.dot=T, draw.l=T, lty.m="dotted", lwd.m=1, invert=F,rectcoords=c(dt[i,2]-dt[i,5],dt[i,2]+dt[i,5]), col=rgb(0,0,0,1))
    if(make.blank)
      labs[i] <- ""
    else
      labs[i] <- dt$categs.cols.names[i]
  }
  
  
  
  axis(1,at=x, labels=labs, las=2)
  
  if(pdf.flag){
    dev.off()
  }
  
}





# EXECUTION ---------

#TEST LINES

# ...
data <-make.data(filepath, categs,sl, categs.names, skip=skip, excluding.fields = excluding.fields, factorings =factorings, delete.outliers=delout, delete.highsd=delhighsd)

CTstats <- make.CTstats(data, alpha=alpha, stcurve.pos=stcurve.pos)
Qstats <- make.Qstats(data, alpha=alpha, stcurve.pos=stcurve.pos)

# save stats
save(data,CTstats,Qstats,file=paste(filepath.res,"_",sep=""))





# GRAPHICS ----------
make.plot(CTstats, pdf.flag = T, paste(filepath.res,"_Ct.pdf",sep="") ,selcolx="categs.cols.names", pad=0.8, wgrap=wgrap, wh=wh, yylab=expression(C[T]))
make.plot(CTstats, pdf.flag = T, paste(filepath.res,"_Ct_blank.pdf",sep="") ,selcolx="categs.cols.names", pad=0.8, wgrap=wgrap, wh=wh, yylab="", make.blank=TRUE)

make.plot(Qstats, pdf.flag = T, paste(filepath.res,"_Q.pdf",sep="") ,selcolx="categs.cols.names", pad=0.8, wgrap=wgrap, wh=wh, yylab="Q")
make.plot(Qstats, pdf.flag = T, paste(filepath.res,"_Q_blank.pdf",sep="") ,selcolx="categs.cols.names", pad=0.8, wgrap=wgrap, wh=wh, yylab="", make.blank=TRUE)


make.plot(CTstats, paste(filepath.res,"_Ct.svg",sep="") ,selcolx="categs.cols.names", pad=0.8, pdf.flag=F, wgrap=wgrap, wh=wh, yylab=expression(C[T]))
make.plot(CTstats, paste(filepath.res,"_Ct_blank.svg",sep="") ,selcolx="categs.cols.names", pad=0.8, pdf.flag=F, wgrap=wgrap, wh=wh, yylab="", make.blank=TRUE)

make.plot(Qstats, paste(filepath.res,"_Q.svg",sep="") ,selcolx="categs.cols.names", pad=0.8, pdf.flag=F, wgrap=wgrap, wh=wh, yylab="Q")
make.plot(Qstats, paste(filepath.res,"_Q_blank.svg",sep="") ,selcolx="categs.cols.names", pad=0.8, pdf.flag=F, wgrap=wgrap, wh=wh, yylab="", make.blank=TRUE)





