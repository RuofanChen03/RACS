# RACS pipeline -- interGENIC regions
# utilities file for interGeneRegions.R script


#########################################################################

reshapeTable <- function(origDATA){

lscaffold <- c()
lbregion <- c()
leregion <- c()

# reshape the table to FIX order...
for (i in origDATA$geneSCFFLD) {
    lscaffold <- c(lscaffold, strsplit(i,':')[[1]][1])
    region <- strsplit(i,':')[[1]][2]
    lbregion <- c(lbregion, as.numeric(strsplit(region,'-')[[1]][1]))
    leregion <- c(leregion, as.numeric(strsplit(region,'-')[[1]][2]))
}
tmpTable <- data.frame(lscaffold,lbregion,leregion)

# the sorting will NOT work for combinations of strings and numbers with an arbitrary number of numnerical digits
#sortedTABLE <- tmpTable[order(tmpTable$lscaffold,tmpTable$lbregion),]
#
# hence we need a more robust and generic approach...
# I) using R-base functions
scfSTRG <- as.numeric(gsub("[^[:digit:]]", "", lscaffold))
names(scfSTRG) <- seq_along(scfSTRG)
sortedTABLE <- tmpTable[as.numeric(names(sort(scfSTRG))),]
# II) alternatively one could use an auxiliary library "stringr"
# eg. str_sort(tmpTable$lscaffold, numeric = TRUE)

# print the shape of the table
cat("DEBUG:: Reshaped table with",dim(sortedTABLE)[1],"rows and",dim(sortedTABLE)[2],"columns",'\n')

return(sortedTABLE)
}
#######################################################################


#########################################################################

dumpData <- function(scaffold,region1,region2, flag) {
   r1 <- as.numeric(region1)+1
   if (region2 != 'xxx') {
       #DBG output: print(paste("R1/R2:::", region1,region2))
       r2 <- as.numeric(region2)-1
       geneSz <- r2-r1+1
      } else {
          r2 <- region2
          geneSz <- 'XXX'
        }

   strng0 <- paste(scaffold,':',r1,'-',r2,sep='')
   strng1 <- paste(strng0,'\t',geneSz,sep='')
   #lstscfld <<- c(lstscfld, strng0)
   #lstreg1 <<- c(lstreg1, r1)
   #lstreg2 <<- c(lstreg2, r2)
   #lstsize <<- c(lstsize, geneSz)
   
   if (flag==1) {
      cat(strng1)
      cat('\n')
   }

   return(c(strng0,r1,r2,geneSz))
}

########################################################################


########################################################################
########################################################################
########################################################################

# OLD intergene regions generating function
interGenesRegionGen_OLD <- function(srtedTABLE) {
SCFLD <- ""
for (i in srtedTABLE$lscaffold) {
    scaffold <- i       #strsplit(i,':')[[1]][1]        #substr(i,1,12)
    region <- strsplit(i,':')[[1]][2]
    print(i)
    if (SCFLD == scaffold) {    #same scaffold
       beginRegion <- strsplit(region,'-')[[1]][1]
       #print(i)
       cat(SCFLD,':',as.numeric(endRegion)+1,'-',as.numeric(beginRegion)-1,'\t',as.numeric(beginRegion)-as.numeric(endRegion))
       cat('\n')
       endRegion <- strsplit(region,'-')[[1]][2]
    } else {    # change scaffold
       if (SCFLD!='') {
           #print(i)
           beginRegion <- strsplit(region,'-')[[1]][1]
           dumpData(SCFLD,endRegion,'xxx')
           dumpData(scaffold,'0',beginRegion)
           endRegion <- strsplit(region,'-')[[1]][2]
       } else { # first scaffold ever... 
              beginRegion <- strsplit(region,'-')[[1]][1]
              endRegion <- strsplit(region,'-')[[1]][2]
              dumpData(scaffold,'0',beginRegion)
         }
      SCFLD <- scaffold
    }
}
# take care of last case...
dumpData(scaffold,endRegion,'xxx')
}

########################################################################
########################################################################
########################################################################


########################################################################





########################################################################

# function to dump the new data generated in a CSV file
saveDATA <- function(data,fileName) {
    write.table(data, file=fileName, sep='\t', row.names=FALSE, quote=FALSE)
}

########################################################################
# Help/Description Fn
errMsgFn <- function(...){

	cat("RACS: intergenic region determination Rscript",'\n')
	cat("----------------------------------------------------",'\n')
	cat("Attention! This script requires 3 arguments!",'\n')
	cat('\t'," i: input file where to read the combined tables from",'\n')
	cat('\t'," ii: reference *gff3* genome file for the organism, eg. 'T_thermophila_June2014.sorted.gff3'", '\n')
	cat('\t'," iii: name of the generated intergenic files, eg. 'interGENs.csv'", '\n')
	cat("---------------------------------------------------",'\n\n')
	cat(paste(...,sep=''),'\n') 
	stop()
}


# Function to process command line arguments
CLArgs <- function(def.ref.file="DATA/T_thermophila_June2014.sorted.gff3",def.out.name='interGENs.csv') {

	checkFile <- function(fileN) {
	   if (file.exists(fileN) == FALSE) { 
		errMsgFn("Error: '",fileN,"' NOT found!")
	   }
	}

# read command-line arguments, expecting three filenames:
#	 i) to read, ii) ref gff3 file for organism, iii) to save results

	args <- commandArgs(trailingOnly=TRUE)

	len.args <- length(args)
	if (len.args < 1) errMsgFn("Error: this scripts requires 3 arguments!")

	# process 1st argument
	inputFile <- args[1]    # eg. "combinedTABLES_BD1-BD2--SORTED"
	checkFile(inputFile)

	if (len.args == 3) {
		refFile <- args[2]	# eg. "DATA/T_thermophila_June2014.sorted.gff3"
		checkFile(refFile)
		outFile <- args[3]
	} else if (len.args == 2 ) {
		   refFile <- args[2]
		   checkFile(refFile)	
		   cat("Assuming default name for output file... '",def.out.name,"'",'\n')
		   outFile <- def.out.name
	       } else {
			errMsgFn("Error: incorrect numnber of arguments!")
		      }

	return(c(inputFile,refFile,outFile))
}


# Function for selecting the input data from a file....
readDATA <- function(filename) {

## read DATA
#combinedDATA=read.csv("~/Downloads/combinedTABLES_BD1-BD2--SORTED",header=F,sep='\t')
#combinedDATA=read.csv("./TABLEE",header=F,sep='\t')
#inputDATA <- read.csv(filename, header=F, sep='\t')

# latest version of RACS-ORF will generate final version of the TABLESwith headers
inputDATA <- read.csv(filename, header=TRUE, sep='')

# give names to the columns...
names(inputDATA)[1]<-"geneSCFFLD"
#names(combinedDATA)[3]<-"Note"
#names(combinedDATA)[4]<-"geneSize"
#names(inputDATA)[5]<-"BD1readsINPUT"
#names(inputDATA)[6]<-"BD1readsIP"
#names(combinedDATA)[7]<-"normBD1readsINPUT"
#names(combinedDATA)[8]<-"normBD1readsIP"
#names(combinedDATA)[9]<-"BD1score"
#names(inputDATA)[10]<-"BD2readsINPUT"
#names(inputDATA)[11]<-"BD2readsIP"
#names(combinedDATA)[12]<-"normBD2readsINPUT"
#names(combinedDATA)[13]<-"normBD2readsIP"
#names(combinedDATA)[14]<-"BD2score"
#names(combinedDATA)[15]<-"BD1BD2enr"

return(inputDATA)
}
#########################################################################
readRefTable <- function(refFile,KWRD='contig') {

library(data.table)

refTableOrig <- read.csv(refFile, header=FALSE, sep='\t')

# exact match
#KWRD='supercontig'	# for T.T.
#refTable  <- data.frame(refTableOrig[refTableOrig$V3==KWRD,]$V1, refTableOrig[refTableOrig$V3==KWRD,]$V5)
# partial match
refTable  <- data.frame(refTableOrig[refTableOrig$V3 %like% KWRD,]$V1, refTableOrig[refTableOrig$V3 %like% KWRD,]$V5)

origSz <- dim(refTable)[1]
# in case that there duplicated
refTable <- unique(refTable)
cat("Original records in ref. table: ", origSz,'\n')
cat("after eliminated possible duplicates...",dim(refTable)[1],'\n')

names(refTable) <- c("scaffold","supercontig")

return(refTable)
}

#########################################################################

scfldCAP <- function(scfld,refTable) {
   return(refTable[refTable$scaffold==scfld,]$supercontig)
}

#########################################################################

warningSCFLD <- function(scfld,begReg,endReg) {
	print(">>>>>>>>>>>>>>>> SUSPICIOUS OVERLYING REGIONs!!!  <<<<<<<<<<<<<<<<")
	cat(scaffold,begReg,endReg,'\n')
	cat(try(lstscfld[length(lstscfld)-1]),try(lstregion1[length(lstscfld)-1]),try(lstregion2[length(lstscfld)-1]),'\n')
	#stop
}

########################################################################




