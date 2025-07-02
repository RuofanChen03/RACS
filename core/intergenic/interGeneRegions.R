#########################################################################
# RACS: InterGenic Regions -- main driver R script
#
# requires "utils_RACS-IGR.R" utilties file with fn defns
#
# part of the det-intergenic.sh pipeline
# used for determing the intergenic regions
#
#
# HOW TO USE this script:
#
# Rscript interGeneRegions.R inputFile refFILE.gff3 outFile
#
#
########################################################################

#
# Function to establish the location of the script... hence the relative path to the utilities file
# this function needs to be in the main driver script!
#

ALLargs <- function() {
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    scrPATHkw <- "--file="
    match <- grep(scrPATHkw, cmdArgs)
    if (length(match) > 0) {
        # Rscript
        scriptLOC <- dirname(normalizePath(sub(scrPATHkw, "", cmdArgs[match])))
    } else {
        # 'source'd via R console
        scriptLOC <- normalizePath(sys.frames()[[1]]$ofile)
    }

    cmdArgs <- commandArgs(trailingOnly =TRUE)
    return(list(scriptLOC,cmdArgs))
}



allArgs <- ALLargs()

print(allArgs)

########################################################################
########################################################################
cat(" RACS v1.0 (2018/2019) -- InterGenic Regions Finder",'\n')
cat("----------------------------------------------------",'\n')
#######################################################################
########################################################################
# load utilities file with fns. defns.

utilitiesFile <- paste(allArgs[[1]],"/utils_RACS-IGR.R",sep='')
print(utilitiesFile)
#print(dirname(sys.frame(1)$ofile)
#sourceDir <- getSrcDirectory(function(dummy) {dummy})
#print(sourceDir)


if (file.exists(utilitiesFile)) {
	source(utilitiesFile)
} else {
	stop("Critical ERRROR: '",utilitiesFile,"' NOT found!")
}

########################################################################

################################


# read files...
files2process <- CLArgs()

inputDATA <- readDATA(files2process[1])


refTable <- readRefTable(files2process[2])


sortedDATA <- reshapeTable(inputDATA)

nbrEntries <- length(sortedDATA$lscaffold)



#for (i in c(1:nbrEntries)) {
#    #a <- sortedDATA[i,][1]	#levels(sortedDATA[i,]$lscaffold)[i]
#    a <-  levels(sortedDATA$lscaffold)[sortedDATA$lscaffold][i]
#    aa <- sortedDATA[i,]$lbregion
#    aaa <- sortedDATA[i,]$leregion
#    cat(a,aa,aaa,'\n')
#}

lstscfld <- c()
lstregion1 <- c()
lstregion2 <- c()
lstSize <- c()
#nbrEntries <- 20

SCFLD <- ""
for (i in c(1:nbrEntries)) {
    scaffold <- levels(sortedDATA$lscaffold)[sortedDATA$lscaffold][i]
    begReg <- sortedDATA[i,]$lbregion
    endReg  <- sortedDATA[i,]$leregion
    cat(i, scaffold,begReg,endReg,'\n')
    if (SCFLD == scaffold) { 	#same scaffold
	#if ((lstscfld[length(lstscfld)]==scaffold)&(begReg<lstregion2[length(lstscfld)])) {
#	 if (begReg>lstregion1[try(length(lstscfld))]) {
#	if (begReg>endReg) {
	if  (begReg > sortedDATA[i-1,]$leregion) {
       beginRegion <- begReg
       #dumpData(SCFLD,endRegion,beginRegion, lstscfld,lstregion1,lstregion2,lstSize, 1)
       myScfld <- dumpData(SCFLD,endRegion,beginRegion, 1)
        lstscfld <- c(lstscfld, myScfld[1])
        lstSize <- c(lstSize, myScfld[4])
        lstregion1 <- c(lstregion1,endRegion+1)
        lstregion2 <- c(lstregion2,beginRegion-1)
       endRegion <- endReg
	} else {	# will to skip this region...
		print(">>>>>>>>>>>>>>>> SUSPICIOUS OVERLYING REGIONs!!!  <<<<<<<<<<<<<<<<")
		cat(scaffold,begReg,endReg,'\n')
		cat(try(lstscfld[length(lstscfld)-1]),try(lstregion1[length(lstscfld)-1]),try(lstregion2[length(lstscfld)-1]),'\n')
		#stop
		}
    } else { 	# change scaffold
      if (SCFLD != "") {
         beginRegion <- begReg
         #dumpData(SCFLD,endRegion,'xxx', lstscfld,lstregion1,lstregion2,lstSize, 1)
         scfCAP <- scfldCAP(SCFLD,refTable)
         #DBG output: print(paste("XXX >~>~>",endRegion,"/",endReg," - ",scfCAP))
         # cautionary guard: delaing with the case of  scfCAP  being empty ""  --> Integer(0)!!!
         #  in principle due to data inconsistencies, keep it for robustness of the code
         if (length(scfCAP)>0) {
		if (endRegion < scfCAP) {	# dealing with cases where the gene ends at the limit of the scaffold...
                   myScfld <- dumpData(SCFLD,endRegion,scfCAP+1, 1)
                   lstscfld <- c(lstscfld, myScfld[1])
                   lstSize <- c(lstSize, myScfld[4])
                   lstregion1 <- c(lstregion1,endRegion+1)
                   lstregion2 <- c(lstregion2,scfCAP)
                } else {
			print(">>>>>>>>>>>>>>>> SUSPICIOUS OVERLYING REGIONs!!!  <<<<<<<<<<<<<<<<")
			cat(scaffold,begReg,endReg,'\n')
			cat(try(lstscfld[length(lstscfld)-1]),try(lstregion1[length(lstscfld)-1]),try(lstregion2[length(lstscfld)-1]),'\n')
			#stop
           	}
                #dumpData(scaffold,'0',beginRegion, lstscfld,lstregion1,lstregion2,lstSize, 1)
                myScfld <- dumpData(scaffold,'0',beginRegion, 1)
                lstscfld <- c(lstscfld, myScfld[1])
                lstSize <- c(lstSize, myScfld[4])
                lstregion1 <- c(lstregion1,0+1)
                lstregion2 <- c(lstregion2,beginRegion-1)
                endRegion <- endReg
	 } else {
		# potential data inconsistency...
		cat("Potential data inconsistency issue detected... please verify your data integrity...",'\n\n')
        }
      } else { 	# first scaffold ever ...
              beginRegion <- begReg	#strsplit(region,'-')[[1]][1]
              endRegion <- endReg	#strsplit(region,'-')[[1]][2]
              #dumpData(scaffold,'0',beginRegion, lstscfld,lstregion1,lstregion2,lstSize, 1)
              myScfld <- dumpData(scaffold,'0',beginRegion, 1)
               lstscfld <- c(lstscfld, myScfld[1])
               lstSize <- c(lstSize, myScfld[4])
               lstregion1 <- c(lstregion1,0+1)
               lstregion2 <- c(lstregion2,beginRegion-1)
         }
      SCFLD <- scaffold
    }
	if (length(lstregion1) != length(lstregion2)) stop()
}
# take care of last case...
#dumpData(scaffold,endRegion,'xxx', lstscfld,lstregion1,lstregion2,lstSize, 1)
 scfCAP <- scfldCAP(scaffold,refTable)
 myScfld <- dumpData(scaffold,endRegion,scfCAP+1, 1)
  lstscfld <- c(lstscfld, myScfld[1])
  lstSize <- c(lstSize, myScfld[4])
 lstregion1 <- c(lstregion1,endRegion+1)
 lstregion2 <- c(lstregion2,scfCAP)

interGenes <- data.frame(lstscfld,lstregion1,lstregion2,lstSize)
names(interGenes) <- c("scaffold",'beggining','end','size')
saveDATA(interGenes,files2process[3])




