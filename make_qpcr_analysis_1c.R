
#
# Regular analysis of qPCR in A. Sanz' lab. 
# Technical replicas are reduced to its mean, their error discarded.
# 
# Biological replicas are then compared using a t-test and a KS-test.
# 
# @AUTHOR: Maximo Sanchez-Aragon
# @EMAIL: maxsanara@gmail.com
#

args = commandArgs(trailingOnly=TRUE) 

# PARAMS -------------
# experimental folders
wd <- args[1] #"/home/maxkavian"
data.fld <- paste(wd,args[2],sep="/") #paste(wd,"qPCR_data",sep="/")               # arg 2
results.fld <- paste(wd,args[3],sep="/") #paste(wd,"qPCR_analysis",sep="/")        # arg 3
data.exp.fld <- paste(c(args[4]),sep="/") #paste(c("20230705"),sep="/")        # arg 4


# experimental filename
data.exp.file <- c(args[5]) 

# housekeeping reference
data.expref.file <- c(args[6]) 

# categories ranks that map categs vector and categs.names vector. These come from the plate reader itself
sl <-  c(args[7]) #  c("a1:a12, b1:d4, b5:d8, e1:g4, e5:g8") arg 7
categs.arg <- args[8] #"st.curve,A1_4,A5_8,D1_4,D5_8" "st.curve,A1_4,A5_8,D1_4,D5_8"   # arg 8
categs.names.arg <- c(args[9]) #c("Standard Curve,tubGS>shRFeSp1 mef(-),tubGS>shRFeSp1 mef(+),tubGS/w[1118] mef(-),tubGS/w[1118] mef(+)")  # arg 9

# OTHER params
stcurve.pos <-  as.numeric(args[10]) # as.numeric("1")  # arg 10 as.numeric(args[10])

skip <- as.numeric(args[11])     # arg 11, lines in th csv that will be skipped by default ***** CURRENTLY DEPRECATED


alpha <- as.numeric(args[12])     # arg 12, alpha for confidence intervals
wgraph <-  as.numeric(args[13])   # arg 13
wh.arg <-  args[14]               # arg 14, width and height of the output graph
delout <-  as.logical(args[15])   # arg 15, should outliers be removed?
delhighsd <- as.logical(args[16]) # arg 16 should high SD points be removed?

alpha <-  as.numeric(args[17]) #arg 17
#yylimCT<- c(-1,2)
#yylimRT<- c(0,2)

# excluding columns from QuatumStudio3 format
excluding.fields.arg <-   args[18]  # "1,3,7,8,10,11,13,14,25,30,31,32,33,29" arg 18, fields excluded from data... 27 HIGHSD # 28 OUTLIERS
# factoring columns in the file (i.e. make them into factors)
factorings.arg <-  args[19]  # "Target.Name,Task,Sample.Name,Well.Position" #  arg 19, columns that need refactoring for R compatibility


# static variables
selector.function <- "mean" # how sample replicas should be represented by?, i.e. centered measurement.





# FUNCTIONS ----------------
ci.bound <- function(x, alpha=.05, low.bound=FALSE){
  n <- length(x)
  a <- alpha/2
  if(low.bound)
    a <- 1-alpha/2
  ci <- sqrt((1/(n-1))*var(x)) * qt(a,n-1)
  
  mean(x) + ci
}
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
  
  if(delete.highsd)
    data <- subset(data,HIGHSD!="Y")
  
  if(delete.outliers)
    if(length(which(names(data) == "OUTLIERRG") ) > 0)
      data <- subset(data,OUTLIERRG!="Y")
  
  n.sample <- aggregate(data, list(data$Sample.Name), length)[,c(1,2)]
  names(n.sample) <- c("Sample.Name","n")
  data <-merge(data,n.sample,by="Sample.Name")
  
  
  
  data
}
keep.lines.f <- function(data,sample.colname = "Sample.Name"){
  lvs <- levels(data[[sample.colname]])[which( nchar(as.character(levels(data[[sample.colname]]))) > 0 )] # excludes empty characters, e.g. standard curve points
  keep.lines <- c()
  for(i in 1:length(lvs))
    keep.lines <- c(keep.lines, min(which(lvs[i]==data[[sample.colname]])))
  data[keep.lines,]
}
sd.f <- function(x){
  sqrt( (length(x)/(length(x)-1))*var(x) )
}
shapiro.pval <- function(x){
  r <- shapiro.test(x)
  r$p.value
}
q.f <- function(x, p=.5){
  quantile(x,probs=p)
}
find.qlevels <- function(plvs=c(.01,.25,.5,.75,.99),tag="Ratio.Q"){
  
  # finds the positions of a tagged list, serving as column names, e.g. quantiles
  
  plvs.n <- paste(tag,plvs*100, sep="")
  lvsi <- c()
  for(i in 1:length(plvs)){
    lvsi[i] <- which(names(data.comb1.stats) == plvs.n[i])
  }
  lvsi
}
make.stats.table <- function(data.comb1,
                             vr ="Ratio",
                             alpha=.05,
                             fncs =c("mean","sd.f","length"),
                             fncs.names=paste(vr,fncs,sep="."),#c("Ratio.Mean","Ratio.SD"),
                             p = c(.01,.05,.25,.5,.75,.95,.99)){
  
  
  
  tbs <- list()
  pn <- as.character(paste(vr,".Q",p*100,sep=""))
  for(i in 1:length(fncs)){ # makes other stats, mean, sd, etc
    tbs[[i]] <- aggregate(subset(data.comb1,select=vr), list(data.comb1$categs.cols), fncs[i])#aggregate(Ratio~categs.cols ,data.comb1, fncs[i])
    names(tbs[[i]])[1] <- "categs.cols"
    names(tbs[[i]])[which(names(tbs[[i]])==vr)] <- fncs.names[i]
  }
  
  # CI interval for this table
  i <- i +1
  tbs[[i]] <- aggregate(subset(data.comb1,select=vr), list(data.comb1$categs.cols),   ci.bound, low.bound=FALSE, alpha=alpha)
  names(tbs[[i]])[1] <- "categs.cols"
  names(tbs[[i]])[which(names(tbs[[i]])==vr)] <-  paste(vr,".CI.low",sep="")
  
  i <- i +1
  tbs[[i]] <- aggregate(subset(data.comb1,select=vr), list(data.comb1$categs.cols),   ci.bound, low.bound=TRUE, alpha=alpha)
  names(tbs[[i]])[1] <- "categs.cols"
  names(tbs[[i]])[which(names(tbs[[i]])==vr)] <- paste(vr,".CI.high",sep="")
  
  
  for(i in 1:length(p)){ # makes the quantiles
    j <- i+length(fncs)+2
    tbs[[j]] <- aggregate(subset(data.comb1,select=vr), list(data.comb1$categs.cols), q.f, p=p[i])#aggregate(Ratio~categs.cols ,data.comb1, q.f, p=p[i])
    names(tbs[[j]])[1] <- "categs.cols"
    names(tbs[[j]])[which(names(tbs[[j]])==vr)] <- pn[i]
  }
  rm(i,j)
  
  R <- merge(tbs[[1]],tbs[[2]], by="categs.cols")
  for(i in 3:length(tbs)){
    R <- merge(R,tbs[[i]], by="categs.cols")
  }
  rm(tbs)
  R
}
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
make.plot <- function(data.comb1.stats, filepath=NULL, wh=c(2,6), vr="Ratio",make.blank=FALSE,  pad =.8, pdf.flag=TRUE, svg.flag=TRUE, yylim=NULL, ...){
  
  
  a1<-paste(vr,".CI.low",sep="")
  a2<-paste(vr,".CI.high",sep="")
  b1<-paste(vr,".sd.f",sep="")
  m1<-paste(vr,".mean",sep="")
  ai1<- which(names(data.comb1.stats) == a1)
  ai2<- which(names(data.comb1.stats) == a2)
  bi1<- which(names(data.comb1.stats) == b1)
  mi1<- which(names(data.comb1.stats) == m1)
  
  #vr <- "Ratio"
  #l<-find.qlevels(plvs=p,tag=paste(vr,".Q",sep=""))
  selcolx <- "categs.cols.names"
  dt <- data.comb1.stats
  
  lvs <- dt[[selcolx]]
  x <- 1:length(lvs)
  
  #yylim <- c(min(dt[,l]), max(dt[,l]))
  if(is.null(yylim))
    yylim <- c(min(dt[,ai1]), max(dt[,ai2]))
  
  
  if(pdf.flag){
    pdf(file=filepath, wh[1],wh[2])
  }else if (svg.flag){
    svg(file=filepath, wh[1],wh[2])
  }
  
  # graphical
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
       xlab="", xaxt="n", ...)
  
  labs<- c()
  for(i in 1:nrow(dt)){
    #draw.point(i,dt[i,l[3]],dt[i,l[1]],dt[i,l[5]],wgraph, draw.dot=T, draw.l=T, lty.m="dotted", lwd.m=1, invert=F,rectcoords=c(dt[i,l[2]],dt[i,l[4]]), col=rgb(0,0,0,1))
    draw.point(i,dt[i,mi1],dt[i,ai1],dt[i,ai2],wgraph, draw.dot=T, draw.l=T, lty.m="dotted", lwd.m=1, invert=F,rectcoords=c(dt[i,mi1]-dt[i,bi1],dt[i,mi1]+dt[i,bi1]), col=rgb(0,0,0,1))
    if(make.blank)
      labs[i] <- ""
    else
      labs[i] <- as.character(dt$categs.cols.names[i])
  }
  
  axis(1,at=x, labels=labs, las=2)
  
  if(pdf.flag || svg.flag){
    dev.off()
  }
  
}



# param processing
wh <- as.numeric(strsplit(wh.arg,",")[[1]]) # width and height for the graphical output
categs <- strsplit(categs.arg,",")[[1]]
categs.names <- strsplit(categs.names.arg,",")[[1]]

basename<- strsplit(as.character(data.exp.file),".", fixed=T)[[1]][1]
basenameref<- strsplit(as.character(data.expref.file),".", fixed=T)[[1]][1]

filepath <- paste(data.fld,data.exp.fld,data.exp.file,sep="/")
filerefpath <- paste(data.fld,data.exp.fld,data.expref.file, sep="/")

results.filename <- paste("GeneExpRes",basename,"_x_",basenameref, sep="")
filepath.res <- paste(results.fld,data.exp.fld,results.filename,sep="/")

excluding.fields <- as.integer(strsplit(excluding.fields.arg,",")[[1]]) # 27 HIGHSD # 28 OUTLIERS
factorings <- as.character(strsplit(factorings.arg,",")[[1]])

# SKIP lines automatic selection 
skip.selector.string <- "Well Position"
skip <- as.integer(system(paste("cat ",filepath," | grep -n \"",skip.selector.string,"\" | cut -d: -f1",sep=""),intern=TRUE)) -1
skip_ref <- as.integer(system(paste("cat ",filerefpath," | grep -n \"",skip.selector.string,"\" | cut -d: -f1",sep=""),intern=TRUE)) -1


data <-make.data(filepath, categs,sl, categs.names, skip=skip, excluding.fields = excluding.fields, factorings =factorings, delete.outliers=delout, delete.highsd=delhighsd)

dataref <-make.data(filerefpath, categs,sl, categs.names, skip=skip_ref, excluding.fields = excluding.fields, factorings =factorings, delete.outliers=delout, delete.highsd=delhighsd)




# builds data table
data1 <- subset(data, Sample.Name!="", select=c(Sample.Name,CT,Quantity,categs.cols, categs.cols.names,n))
dataref1 <- subset(dataref, Sample.Name!="", select=c(Sample.Name,CT,Quantity,n))

rm(data, dataref)

names(dataref1)[which(names(dataref1) == "Quantity")] <- "QuantityRef"
names(dataref1)[which(names(dataref1) == "CT")] <- "CTRef"
names(dataref1)[which(names(dataref1) == "n")] <- "nRef"

# makes the means, medians, etc of replicas
dt1 <- aggregate(cbind(CT,Quantity)~Sample.Name, data1, selector.function)
data1 <- merge(dt1, subset(keep.lines.f(data1),select=c(Sample.Name,categs.cols,categs.cols.names)))
rm(dt1)

dtref1 <- aggregate(cbind(CTRef,QuantityRef)~Sample.Name, dataref1, selector.function)
dataref1 <- merge(dtref1, subset(keep.lines.f(dataref1),select=c(Sample.Name)))
rm(dtref1)


# combines every possible item from Sample.Name in data1 to every possible item from Sample.Name to dataref1
# NOTE: theoretically only one biological sample should remain in the data table.
data.comb <- merge(data1,dataref1, by="Sample.Name")

rm(data1,dataref1)

data.comb1 <- data.frame(data.comb, Ratio=data.comb$Quantity/data.comb$QuantityRef, DeltaCT=data.comb$CT - data.comb$CTRef)
rm(data.comb)

data.comb1$Sample.Name <- factor(data.comb1$Sample.Name)
data.comb1$categs.cols <- factor(data.comb1$categs.cols)






# MAKES the stats table ---------------------------
data.comb1.stats <- merge(make.stats.table(data.comb1, vr="Ratio", alpha=alpha), make.stats.table(data.comb1, vr="DeltaCT", alpha=alpha), by="categs.cols")
pretb <- keep.lines.f(subset(data.comb1, select=c(categs.cols,categs.cols.names) ), sample.colname = "categs.cols")

data.comb1.stats <- merge(pretb, data.comb1.stats)
rm(pretb)




# GRAPHICAL output -----------------------------
make.plot(data.comb1.stats, filepath=paste(filepath.res,"_Ratio.pdf",sep=""), vr="Ratio", ylab="Ratio", yylim=NULL)
make.plot(data.comb1.stats, filepath=paste(filepath.res,"_Ratio.svg",sep=""), vr="Ratio", ylab="Ratio", yylim=NULL, pdf.flag = FALSE)

make.plot(data.comb1.stats, filepath=paste(filepath.res,"_DeltaCT.pdf",sep=""), vr="DeltaCT", yylim=NULL, ylab=expression(paste(Delta,C[t],sep="")))
make.plot(data.comb1.stats, filepath=paste(filepath.res,"_DeltaCT.svg",sep=""), vr="DeltaCT", yylim=NULL, ylab=expression(paste(Delta,C[t],sep="")), pdf.flag = FALSE)


# STATISTICAL comparisons --------------------------
make.ks.table <- function(data.comb1,comparer="Ratio", testf="ks"){
  data.comb1.ls <- split(data.comb1, f=data.comb1$categs.cols)
  
  p.vals <- c()
  D.vals <- c()
  categs.vals <- c()
  categsC.vals <- c()
  categs.vals.names <- c()
  categsC.vals.names <- c()
  for(i in 1:(length(data.comb1.ls)-1)) {
    for(j in (i+1):(length(data.comb1.ls)) ){
  
      # Test line
      if(testf == "ks")
        R <-ks.test(data.comb1.ls[[i]][[comparer]], data.comb1.ls[[j]][[comparer]])
      else
        R <-t.test(data.comb1.ls[[i]][[comparer]], data.comb1.ls[[j]][[comparer]])
      
      p.vals <- c(p.vals, R$p.value)
      D.vals <- c(D.vals, R$statistic)
      categs.vals <- c(categs.vals, names(data.comb1.ls)[i])
      categsC.vals <- c(categsC.vals,names(data.comb1.ls)[j])
      categs.vals.names <- c(categs.vals.names, as.character(data.comb1.ls[[i]]$categs.cols.names[1]))
      categsC.vals.names <- c(categsC.vals.names, as.character(data.comb1.ls[[j]]$categs.cols.names[1]))
    }
  }
  
  results <- data.frame(categs.vals,categsC.vals,categs.vals.names, categsC.vals.names, p.vals,D.vals)
  attr(results, "description") <- "KS tests on several categories"
  
  results
}
ks.table.ratio <- make.ks.table(data.comb1, comparer="Ratio")
ks.table.deltaCT <- make.ks.table(data.comb1, comparer="DeltaCT")
t.table.ratio <- make.ks.table(data.comb1, comparer="Ratio", testf="t.test")
t.table.deltaCT <- make.ks.table(data.comb1, comparer="DeltaCT", testf="t.test")


# SAVE
attr(data.comb1,"alpha") <- alpha
attr(data.comb1.stats,"alpha") <- alpha
save(data.comb1,data.comb1.stats, 
     ks.table.deltaCT, ks.table.ratio,
     t.table.deltaCT, t.table.ratio, 
     file=paste(filepath.res,"_b",sep=""))





# LATEX printing ------------------------
#wd <- "/home/maxkavian"
#results.fld <- paste(wd,"qPCR_analysis",sep="/")        
#data.exp.fld <- c("20230705")

efile <- paste(results.filename,"_b",sep=""); #"GeneExpResMax_RFeSp1_2023-07-09-080800_Results_x_Max_Act88f_rep1_2023-07-05_080147_Results_b"
path <- paste(results.fld,data.exp.fld,efile,sep="/")

setwd(results.fld)
source("qpcr_funcs.R")


# LOAD --------------------
load(path)

# FUNCTIONS for latex table in data.comb1.stats -----------------
find.names.indexes <- function(data.comb1.stats, selname="Ratio"){
  pref <- c()
  for(i in 1:length(names(data.comb1.stats)))
    pref <- c(pref, strsplit(names(data.comb1.stats),".", fixed=TRUE)[[i]][1])
  
  indexes<-which(pref == selname)
  indexes
}
make.pretable <- function(data.comb1.stats,selname="Ratio", excluded.cols=c("length")){
  
  # data.comb1.stats : data table
  # selname          : selector calculation
  # excluded.cols    : excluded columns (suffixes)
  
  indexes <- find.names.indexes(data.comb1.stats,selname=selname)
  dt <- subset(data.comb1.stats, select=c("categs.cols.names",names(data.comb1.stats)[indexes]) )
  
  for(ex in excluded.cols)
    dt <- dt[,-which(names(dt)==paste(selname,".",ex,sep="") )]
  
  dt
}
make.latex.column <- function(data, n, resolution=3, ...){
  
  cl <- strsplit(as.character(data[,n]), ".", fixed=TRUE, ...)
  cl1 <- c()
  cl2 <- c()
  for(i in 1:length(cl)){
    cl1 <- c(cl1, cl[[i]][1])
    cl2 <- c(cl2, substring(cl[[i]][2],0,resolution) )
  }
  
  R <- data.frame(cl1,cl2)
  names(R) <- c(paste(names(data)[n],"_p",sep=""), paste(names(data)[n],"_s",sep=""))
  R
}

# EXECUTION ------------
excluded.cols <- c("length")

print.latex.table(make.pretable(data.comb1.stats, selname="DeltaCT", excluded.cols=c("length")), path, excluded.cols=c(), text.cols = "categs.cols.names", suffix="_DATACOMB_STATS_DeltaCT", dig=2)

print.latex.table(make.pretable(data.comb1.stats, selname="Ratio", excluded.cols=c("length")), path, excluded.cols=c(), text.cols = "categs.cols.names", suffix="_DATACOMB_STATS_Ratio", dig=2)

print.latex.table(data.comb1, path, excluded.cols=c("Sample.Name","categs.cols"), text.cols=c("categs.cols.names"), suffix="_DATACOMB", dig=2)

print.latex.table(ks.table.deltaCT, path, excluded.cols=c("categs.vals","categsC.vals"), text.cols = c("categs.vals.names","categsC.vals.names"), suffix="_KS_TEST_DeltaCT", dig=2)
print.latex.table(ks.table.ratio, path, excluded.cols=c("categs.vals","categsC.vals"), text.cols = c("categs.vals.names","categsC.vals.names"), suffix="_KS_TEST_Ratio", dig=2)
print.latex.table(t.table.deltaCT, path, excluded.cols=c("categs.vals","categsC.vals"), text.cols = c("categs.vals.names","categsC.vals.names"), suffix="_T_TEST_DeltaCT", dig=2)
print.latex.table(ks.table.ratio, path, excluded.cols=c("categs.vals","categsC.vals"), text.cols = c("categs.vals.names","categsC.vals.names"), suffix="_T_TEST_Ratio", dig=2)















