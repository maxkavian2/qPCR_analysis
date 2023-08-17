
#
# This file processes the data from a specific experiment to find the fold increase of expression
# of a gene relative to a housekeeping (constant)-expressed gene.
# 
# The script finds the distribution of the means of all per category replicas (technical and biological)
# through a t-Student distribution, and performs a Montecarlo sampling of the differences to compute
# the DeltaDeltaCT parameter, using the gene expression reference provided by the user. The resulting 
# distribution of relative fold increases are subjected to KS-analysis afterwards across all possible pairs.
#
# It requires the execution of make_qpcr_analysis1.R first.
# 
# @AUTHOR: Maximo Sanchez-Aragon
# @EMAIL: maxsanara@gmail.com
# 

args = commandArgs(trailingOnly=TRUE) 

# PARAMS -------
wd <- args[1] #"/home/maxkavian/MAX_DOCUMENTS_debian/UOG_project_sep2022" # arg 1
data.fld <- paste(wd,args[2],sep="/")   #paste(wd,"qPCR_data",sep="/")           # arg 2
results.fld <- paste(wd,args[3],sep="/") #paste(wd,"qPCR_analysis",sep="/")    # arg 3
data.exp.fld <- args[4]   #c("20230705")                       # arg 4

eg.filename <- args[5] #"Max_RFeSp1_2023-07-09-080800_Results_" # tested gene, arg 5
rg.filename <- args[6] #"Max_Act88f_rep1_2023-07-05_080147_Results_" # reference gene results, arg 6


exp.lev.arg <- args[7] #c("A1_4,A5_8,D5_8,A5_8") # one-on-one mapping of experimental and control, arg 7
ctr.lev.arg <- args[8] #c("D1_4,D5_8,D1_4,A1_4") # arg 8

monte.carlo.res.arg <- args[9] #"1e+4" # number of items of the monte carlo simulation, arg 9

prob.vals.arg <- args[10] #".01,.05,.25,.5,.75,.95,.99" # arg 10
wh.arg<- args[11] #"2,8" # arg 11
wgrap.arg<- args[12] #"0.8" # arg 12

# KS tests 
lve.arg <- args[13] #"A1_4xD1_4,A1_4xD1_4,A1_4xD1_4,A5_8xD5_8,A5_8xD5_8,A5_8xA1_4" # arg 13
lvc.arg <- args[14] #"A5_8xD5_8,D5_8xD1_4,A5_8xA1_4,D5_8xD1_4,A5_8xA1_4,D5_8xD1_4" # arg 14

dig <- as.integer(args[15]) # number of digits




# Pre-processed arguments
lve <- strsplit(lve.arg,",")[[1]]
lvc <- strsplit(lvc.arg,",")[[1]]
monte.carlo.res <- as.numeric(monte.carlo.res.arg)
results.filename <- paste("GeneExpRes",eg.filename,"x_",rg.filename, sep="")
prob.vals <- as.numeric(strsplit(prob.vals.arg,",")[[1]]) # values from the bootstrap used to get the distribution.
exp.lev <- as.character(strsplit(exp.lev.arg,",")[[1]]) # one-on-one mapping of experimental and control
ctr.lev <- as.character(strsplit(ctr.lev.arg,",")[[1]])

wh<-as.numeric(strsplit(wh.arg,",")[[1]])
wgrap<-as.numeric(wgrap.arg)






# sanity checks 
if(length(exp.lev) != length(ctr.lev))
  stop["Max says: experimental and control categories should be of the same length."]

rg.path <- paste(results.fld,data.exp.fld,rg.filename, sep="/")
eg.path <- paste(results.fld,data.exp.fld,eg.filename, sep="/")

if(!file.exists(rg.path))
  stop("[Max] says: the reference gene file is missing, exiting ...")

if(!file.exists(eg.path))
  stop("[Max] says: the experimental gene file is missing, exiting ...")



# FUNCTIONS ----------------------------
make.montecarlo.t.diff <- function(n, dgfs, means, sems){
  
  # It returns a sample of n elements, which are the random difference 
  # between two t-student distributed variables, with degrees of freedom dgfs[1]
  # and dgfs[2], with means means[1] and means[2], with standard errors of the mean
  # sems[1] and sems[2]

  r <- means[2]-rt(n, dgfs[2])*sems[2] - (means[1]-rt(n, dgfs[1])*sems[1])
  r
}
make.montecarlo.diff <- function(v1,v2, p=c(.01,.05,.25,.5,.75,.95,.99), return.sample=FALSE){
  # makes the difference of two vectors of the same size
  # and then computes the specified quantiles.
  s <- sample(v2,length(v2),replace=FALSE)- sample(v1,length(v2),replace=FALSE)
  if(return.sample){
    return(s)}
    
  r <- quantile(s, probs=p)
  
  attr(r, "mean") <- mean(s)
  attr(r, "var") <- var(s)
  attr(r, "sd") <- sqrt(var(s)*length(s)/(length(s)-1))
  
  r
}
make.expression.diff <- function(egCTstats, rgCTstats,...){
  
  #
  # Computes the expression increase using the CT values
  # egCTstats : statistical table for the target gene (e.g. RFeSp)
  # rgCTstats : statistical table for the reference gene (e.g. actin)
  # ... : arguments passed to make.montecarlo.diff, e.g. vector of probabilities.
  #
  
  explevs <- c()
  ctrlevs <- c()
  explevs.nm <- c()
  ctrlevs.nm <- c()
  rs <- c()
  means <- c()
  vars <- c()
  sds <- c()
  
  for(i in 1:length(exp.lev)){
    
    b <- subset(egCTstats, categs.cols==exp.lev[i], select=c(CT.mean,CT.SEM,CT.n,categs.cols))
    a <- subset(rgCTstats, categs.cols==exp.lev[i], select=c(CT.mean,CT.SEM,CT.n,categs.cols))
    
    ba1 <- make.montecarlo.t.diff(monte.carlo.res, 
                                  c(a$CT.n[1]-1,b$CT.n[1]-1),
                                  c(a$CT.mean[1],b$CT.mean[1]),  
                                  c(a$CT.SEM[1],b$CT.SEM[1]))
    
    b <- subset(egCTstats, categs.cols==ctr.lev[i], select=c(CT.mean,CT.SEM,CT.n,categs.cols))
    a <- subset(rgCTstats, categs.cols==ctr.lev[i], select=c(CT.mean,CT.SEM,CT.n,categs.cols))
    
    ba2 <- make.montecarlo.t.diff(monte.carlo.res, 
                                  c(a$CT.n[1]-1,b$CT.n[1]-1),
                                  c(a$CT.mean[1],b$CT.mean[1]),  
                                  c(a$CT.SEM[1],b$CT.SEM[1]))
    rm(a,b)
    
    explevs <- c(explevs, exp.lev[i])
    ctrlevs <- c(ctrlevs, ctr.lev[i])
    explevs.nm <- c(explevs.nm, names(exp.lev)[i])
    ctrlevs.nm <- c(ctrlevs.nm, names(ctr.lev)[i])
    
    r <- make.montecarlo.diff(ba2,ba1, ...)
    
    means <- c(means, attr(r, "mean"))
    vars <- c(vars, attr(r, "var"))
    sds <- c(sds, attr(r, "sd"))
    
    if(i == 1)
      rs <- r
    else
      rs <- rbind(rs,r)
    
  }
  
  R <- data.frame(rs,explevs,ctrlevs,explevs.nm,ctrlevs.nm, means, vars, sds)
  
  row.names(R) <- NULL
  
  R
}
make.random.diff.sample <- function(egCTstats, rgCTstats,...){

  explevs <- c()
  ctrlevs <- c()
  explevs.nm <- c()
  ctrlevs.nm <- c()
  rs <- c()

  for(i in 1:length(exp.lev)){
    
    b <- subset(egCTstats, categs.cols==exp.lev[i], select=c(CT.mean,CT.SEM,CT.n,categs.cols))
    a <- subset(rgCTstats, categs.cols==exp.lev[i], select=c(CT.mean,CT.SEM,CT.n,categs.cols))
    
    ba1 <- make.montecarlo.t.diff(monte.carlo.res, 
                                  c(a$CT.n[1]-1,b$CT.n[1]-1),
                                  c(a$CT.mean[1],b$CT.mean[1]),  
                                  c(a$CT.SEM[1],b$CT.SEM[1]))
    
    b <- subset(egCTstats, categs.cols==ctr.lev[i], select=c(CT.mean,CT.SEM,CT.n,categs.cols))
    a <- subset(rgCTstats, categs.cols==ctr.lev[i], select=c(CT.mean,CT.SEM,CT.n,categs.cols))
    
    ba2 <- make.montecarlo.t.diff(monte.carlo.res, 
                                  c(a$CT.n[1]-1,b$CT.n[1]-1),
                                  c(a$CT.mean[1],b$CT.mean[1]),  
                                  c(a$CT.SEM[1],b$CT.SEM[1]))
    rm(a,b)
    
    
    r <- make.montecarlo.diff(ba2,ba1, return.sample=TRUE, ...)
    
    explevs <- c(explevs, rep(exp.lev[i],length(r)) )
    ctrlevs <- c(ctrlevs, rep(ctr.lev[i],length(r)) )
    explevs.nm <- c(explevs.nm, rep(names(exp.lev)[i],length(r)) )
    ctrlevs.nm <- c(ctrlevs.nm, rep(names(ctr.lev)[i],length(r)) )
    
    
    rs <- c(rs,r)
    
  }
  
  R <- data.frame(rs,explevs,ctrlevs,explevs.nm,ctrlevs.nm)
  
  R
}
make.ks.tests <- function(egCTstats, rgCTstats, lve=NULL,lvc=NULL){
  data <- make.random.diff.sample(egCTstats,rgCTstats)
  com.categ <- paste(data$explevs,"x",data$ctrlevs,sep="")
  com.categ.nm <- paste(data$explevs.nm, " C:",data$ctrlevs.nm,sep="")
  data <- data.frame(data, com.categ,com.categ.nm)
  pv <- c() # p.value
  D <- c()  # D statistic in the KS test
  
  lve.nm <- c()
  lvc.nm <- c()
  
  for(i in 1:length(lve)){
    
    a<- subset(data, com.categ==lve[i])
    b<- subset(data, com.categ==lvc[i])
    
    #make.ks.test(b$rs,y=a$rs, alternative = "less")
    pv <-c(pv, ks.test(2^(-a$rs),2^(-b$rs), alternative="greater")$p.value) #"greater"
    D <-c(D, ks.test(2^(-a$rs),2^(-b$rs), alternative="greater")$statistic)   #"greater"
    lve.nm <- c(lve.nm, data$com.categ.nm[which(data$com.categ == lve[i])[1]])
    lvc.nm <- c(lvc.nm, data$com.categ.nm[which(data$com.categ == lvc[i])[1]])
  }
  
  Rcomp <- data.frame(pv,D,lve,lvc,lve.nm,lvc.nm)
  
  Rcomp
}
make.plot <- function(R, filepath="", wgrap=1, wh=c(2,6), pdf.flag=T, pad=.8, yylab="", plevs= c(0.05,.25,.5,.75,0.95), make.blank=FALSE ,...){
  
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
  par(mar=c(24,6-mar.red,2.5-mar.red,2.5-mar.red)+0.1)
  #par(fig = c(0,1,0,1))
  # sets the size of the character
  par(cex=.9)
  
  csel <- paste("X",as.character(plevs*100),".",sep="")
  
  yylim <- c(min(2^(-R[,csel]),na.rm=TRUE),max(2^(-R[,csel]),na.rm=TRUE))
  
  plot(NULL,type="n",ylim=yylim, xlim=c(1-pad,nrow(R)+pad), xaxt="n",xlab="",...)
  
  for(i in 1:nrow(R)){
    draw.point( i, 2^(-R[[csel[3]]][i]), 2^(-R[[csel[1]]][i]), 2^(-R[[csel[5]]][i]),w=wgrap,rectcoords=c(2^(-R[[csel[2]]][i]),2^(-R[[csel[4]]][i])), col="black" )
  }
  
  labs<-rep("",times=nrow(R))
  if(!make.blank)
    labs=paste(R$explevs.nm," ref:",R$ctrlevs.nm,sep=" ")
  
  axis(1,at=c(1:nrow(R)), labels=labs, las=2)
  
  if(pdf.flag){
    dev.off()
  }
  
}  
make.p.plot <- function(Rcomp, filepath="", wgrap=1, wh=c(2,6), pdf.flag=T, pad=.8, ...){
  
  if(pdf.flag){
    pdf(file=filepath, wh[1],wh[2])
  }else{
    svg(file=filepath, wh[1],wh[2])
  }
  
  w <-wgrap
  vals <- Rcomp$pv
  yylim<-c(min(vals,na.rm=TRUE), max(vals,na.rm=TRUE))
  if(max(yylim) > 1) yylim = c(yylim[1],1)
  
  # gaphical
  mar.red <- 1
  # changes the likelihood to switch to scientific notation
  options(scipen = 5)
  #set the distances of the axes labels to the plot
  par(mgp=c(2.4,1,.4))
  # sets the margins
  par(mar=c(24,6-mar.red,2.5-mar.red,2.5-mar.red)+0.1)
  #par(fig = c(0,1,0,1))
  # sets the size of the character
  par(cex=.9)
  plot(NULL, type="n", ylim=yylim,xlim=c(1-pad,nrow(Rcomp)+pad),xlab="",xaxt="n",... )
  
  for(i in 1:nrow(Rcomp)){
    rect(i-w/2, 0, i+w/2, vals[i], border=NA, col="steelblue4")
    if(vals[i]<1e-2)
      rect(i-w/2, 0, i+w/2, 0.005, border=NA, col="steelblue4")
  }
  
  lines(c(-1,nrow(Rcomp)+1), c(0.01,0.01), lty="dotted", lwd=1)
  lines(c(-1,nrow(Rcomp)+1), c(0.05,0.05), lty="dotted", lwd=1)
  lines(c(-1,nrow(Rcomp)+1), c(0.10,0.10), lty="dotted", lwd=1)
  
  labs=rep("",times=nrow(Rcomp))#paste(Rcomp$lve.nm, "vs.",Rcomp$lvc.nm,sep=" ")
  axis(1,at=c(1:nrow(Rcomp)), labels=labs, las=2)
  
  if(pdf.flag){
    dev.off()
  }
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









# EXECUTION ------------------------------
load(rg.path)
rgCTstats <- CTstats 
load(eg.path)
egCTstats <- CTstats 

rm(CTstats,Qstats,data)

names(ctr.lev) <- as.character(egCTstats$categs.cols.names)[match(ctr.lev, as.character(egCTstats$categs.cols)) ]
names(exp.lev) <- as.character(egCTstats$categs.cols.names)[match(exp.lev, as.character(egCTstats$categs.cols)) ]


R <- make.expression.diff(egCTstats, rgCTstats, p=prob.vals)

attr(R, "monte.carlo.resolution") <- monte.carlo.res

# KS tests on bootstrap samples.
Rcomp <- make.ks.tests(egCTstats,rgCTstats,lve = lve,lvc=lvc)

save(R,Rcomp,file=paste(results.fld,data.exp.fld,results.filename, sep="/"))








# GRAPHICS ------------------
custylab <- expression(paste(2^-Delta,""^Delta,""^C[t] ,sep=""))

make.plot(R, pdf.flag=T, filepath=paste(results.fld,"/",data.exp.fld,"/",results.filename,"_comp.pdf",sep=""), wgrap=wgrap, wh=wh, ylab=custylab)
make.plot(R, pdf.flag=T, filepath=paste(results.fld,"/",data.exp.fld,"/",results.filename,"_comp_blank.pdf",sep=""), wgrap=wgrap, wh=wh, ylab=custylab, make.blank = TRUE)

make.plot(R, pdf.flag=F, filepath=paste(results.fld,"/",data.exp.fld,"/",results.filename,"_comp.svg",sep=""), wgrap=wgrap, wh=wh, ylab=custylab)
make.plot(R, pdf.flag=F, filepath=paste(results.fld,"/",data.exp.fld,"/",results.filename,"_comp_blank.svg",sep=""), wgrap=wgrap, wh=wh, ylab=custylab, make.blank = TRUE)


# KS tests (CDF of x > CDF of y)
make.p.plot(Rcomp, filepath=paste(results.fld,"/",data.exp.fld,"/",results.filename,"_comp_p.pdf",sep=""), pdf.flag=T, wgrap=wgrap, wh=wh, ylab="")

make.p.plot(Rcomp, filepath=paste(results.fld,"/",data.exp.fld,"/",results.filename,"_comp_p.svg",sep=""), pdf.flag=F, wgrap=wgrap, wh=wh, ylab="")



# LATEX tables export -------------------
#wd <- "/home/maxkavian"
#results.fld <- paste(wd,"qPCR_analysis",sep="/")        
#data.exp.fld <- c("20230705")

efile <-results.filename; #"GeneExpResMax_RFeSp1_2023-07-09-080800_Results_x_Max_Act88f_rep1_2023-07-05_080147_Results_"

path <- paste(results.fld,data.exp.fld,efile,sep="/")


setwd(results.fld)
source("qpcr_funcs.R")

# LOAD --------------------
load(path)

print.latex.table(R, path, excluded.cols=c("explevs","ctrlevs"), text.cols=c("explevs.nm","ctrlevs.nm"), suffix="MONTECARLO_DATA_", dig=dig)
print.latex.table(Rcomp, path, excluded.cols=c("lve","lvc"), text.cols=c("lve.nm","lvc.nm"), suffix="MONTECARLO_TESTS_", dig=dig)











