

require(xtable)

print.latex.table <- function(R, filepath, excluded.cols=c("explevs","ctrlevs"), text.cols=c("explevs.nm","ctrlevs.nm"), suffix="_", 
                              dig=3){
  
  # CREATES a generic pre-formatted latex compliant table.
  # the code excludes non-relevant columns, shifts text columns to the beginning and 
  # generate the right table header.
  # 
  # ARGUMENTS 
  # R             : The data table from which values will be selected
  # filepath      : Path to the original data file, used here for results naming purposes.
  # excluded.cols : Columns that will be excluded
  # text.cols     : Text columns that will be shifted to the beginning of the table
  # suffix        : Added text for later identification
  # dig           : number of digits
  
  # It eliminates the excluded columns from the data
  pretab <- R
  if(length(excluded.cols) >= 1)
  for(i in 1:length(excluded.cols))
    pretab <- pretab[,-which(names(pretab) == excluded.cols[i])]
  
  # It stores the position of the text columns and move it to the beginning
  if(length(text.cols) >= 1)
  for(i in 1:length(text.cols)){
    cpos <- which (names(pretab) == text.cols[length(text.cols)-i+1])[1]
    pretab <- pretab[,c(cpos,which(c(1:ncol(pretab)) != cpos)) ] # rearranges the table for latex formatting
  }
  
  align.arg <- c(rep("c|", times=length(text.cols)+1),rep("r@{.}l|", times=(ncol(pretab) - length(text.cols)) ) )
  display.arg <- c(rep("s", times=length(text.cols)+1),rep("f",times=(ncol(pretab) - length(text.cols))) ) 
  
  latex.code <- xtable(pretab, caption="", label="", digits=dig, align=align.arg, display=display.arg )
  
  print(latex.code, file=paste(filepath,suffix,sep="") )
}
