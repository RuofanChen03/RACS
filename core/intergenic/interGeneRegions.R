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

SCFLD <- ""

for (i in 1:nbrEntries) {
  scaffold <- levels(sortedDATA$lscaffold)[sortedDATA$lscaffold][i]
  begReg <- sortedDATA[i, ]$lbregion
  endReg <- sortedDATA[i, ]$leregion
  
  cat(i, scaffold, begReg, endReg, '\n')
  
  if (SCFLD == scaffold) {
    if (begReg > sortedDATA[i - 1, ]$leregion) {
      beginRegion <- begReg
      myScfld <- dumpData(SCFLD, endRegion, beginRegion, 1)
      
      lstscfld <- c(lstscfld, myScfld[1])
      lstSize <- c(lstSize, myScfld[4])
      lstregion1 <- c(lstregion1, endRegion + 1)
      lstregion2 <- c(lstregion2, beginRegion - 1)
      
      endRegion <- endReg
    } else {
      print(">>>>>>>>>>>>>>>> SUSPICIOUS OVERLYING REGIONs!!! <<<<<<<<<<<<<<<<")
      cat(scaffold, begReg, endReg, '\n')
      cat(
        try(lstscfld[length(lstscfld) - 1]),
        try(lstregion1[length(lstscfld) - 1]),
        try(lstregion2[length(lstscfld) - 1]), '\n'
      )
    }
  } else {
    if (SCFLD != "") {
      beginRegion <- begReg
      scfCAP <- scfldCAP(SCFLD, refTable)
      
      if (length(scfCAP) > 0) {
        if (endRegion < scfCAP) {
          myScfld <- dumpData(SCFLD, endRegion, scfCAP + 1, 1)
          lstscfld <- c(lstscfld, myScfld[1])
          lstSize <- c(lstSize, myScfld[4])
          lstregion1 <- c(lstregion1, endRegion + 1)
          lstregion2 <- c(lstregion2, scfCAP)
        } else {
          print(">>>>>>>>>>>>>>>> SUSPICIOUS OVERLYING REGIONs!!! <<<<<<<<<<<<<<<<")
          cat(scaffold, begReg, endReg, '\n')
          cat(
            try(lstscfld[length(lstscfld) - 1]),
            try(lstregion1[length(lstscfld) - 1]),
            try(lstregion2[length(lstscfld) - 1]), '\n'
          )
        }
        
        myScfld <- dumpData(scaffold, '0', beginRegion, 1)
        lstscfld <- c(lstscfld, myScfld[1])
        lstSize <- c(lstSize, myScfld[4])
        lstregion1 <- c(lstregion1, 1)
        lstregion2 <- c(lstregion2, beginRegion - 1)
        endRegion <- endReg
      } else {
        cat("Potential data inconsistency issue detected... please verify your data integrity...", '\n\n')
      }
    } else {
      beginRegion <- begReg
      endRegion <- endReg
      myScfld <- dumpData(scaffold, '0', beginRegion, 1)
      lstscfld <- c(lstscfld, myScfld[1])
      lstSize <- c(lstSize, myScfld[4])
      lstregion1 <- c(lstregion1, 1)
      lstregion2 <- c(lstregion2, beginRegion - 1)
    }
    SCFLD <- scaffold
  }

  if (length(lstregion1) != length(lstregion2)) stop()
}

scfCAP <- scfldCAP(scaffold, refTable)
myScfld <- dumpData(scaffold, endRegion, scfCAP + 1, 1)

lstscfld <- c(lstscfld, myScfld[1])
lstSize <- c(lstSize, myScfld[4])
lstregion1 <- c(lstregion1, endRegion + 1)
lstregion2 <- c(lstregion2, scfCAP)

interGenes <- data.frame(lstscfld,lstregion1,lstregion2,lstSize)
names(interGenes) <- c("scaffold",'beggining','end','size')
saveDATA(interGenes,files2process[3])




