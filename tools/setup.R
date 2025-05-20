#!/usr/bin/Rscript

###
# Script and Functions for seting up required R-packages for RACS comparison tools 
###

# define dependencies, ie. needed to be installed for the code to run
# first packages that can be installaed from CRAN using install.packages()


if (interactive()) {
	# interactive mode
	# check that variable "pckges" is defined
	pckgesVar <- "pckges" %in% ls()
	if (!pckgesVar) {
		message("You need to define 'pckges' to install!")
		message('E.g. pckges <- "plotly"')
		stop("")
	}
} else {
	# default value for packages when running in batch mode
	pckges <- c("xlsx","plotly")

	# allow to define packages to check via CLA
	CLAs <- commandArgs(trailing=TRUE)
	if (length(CLAs) > 0) {
		pckges <- c(CLAs)
	}
}


##############################################################################################
# functions for checking whether dependencies are installed and install all needed packages


# function to install neeeded packages
NeededPackages <- function(pckges, otherPckgs="", def.mirror='https://cloud.r-project.org') {
	RverM <- as.numeric(R.Version()['major'])
	Rverm <- as.numeric(R.Version()['minor'])

        availablePckges <- .packages(all.available = TRUE)
	#print(availablePckges)

	# deal with packages from CRAN
        needTOinstall <- !(pckges %in% availablePckges)
	if (sum(needTOinstall) != 0) {
	    cat("Requested packages from CRAN:")
	    print(pckges)
	    cat("installing...", pckges[needTOinstall], '\n')
	    for (pck in pckges[needTOinstall]) {
                install.packages(pck,repos=def.mirror)
	    }
	}

	# deal with packages from BioConductor
        needTOinstall <- !(otherPckgs %in% availablePckges)
	if ((otherPckgs!="") && (sum(needTOinstall) != 0) ) {
	    cat("Requested packages:")
	    print(otherPckgs)
	    cat("installing from BioConductor...", otherPckgs[needTOinstall], '\n')
            for (pck in otherPckgs[needTOinstall]) {
		print(RverM); print(Rverm)
		if ((RverM >= 3 ) && (Rverm > 5)) {
			# newer R version
                        install.packages("BiocManager", repos=def.mirror)
                        BiocManager::install(pck)
		} else {
			# older R versions...
			source("https://bioconductor.org/biocLite.R") 
			library(BiocInstaller)
			BiocInstaller::biocLite(pck)
		}
            }
	}
}

# function to check the versions of a given set of packages...
checkVersion <- function(pckges) {

        print(sessionInfo())

	print(pckges)

        for (pck in pckges){
                cat(pck, as.character(packageVersion(pck)), '\n')
        }

}

##############################################################################################



# check and install required packages
NeededPackages(pckges)

# display versions of the installed/required packages
checkVersion(pckges)

