#!/bin/bash

# Tool for downloading "Oxytrichia" data from NCBI_SRA
# This script is part of
# RACS v1.0 (2018/2019) -- Open source tools for Analizing ChIP-Seq data                                   

# It requires the "NCBI SRA (Sequence Read Archive)" to download the data
# whcioh can be obtained from	https://github.com/ncbi/sra-tools

# if it is not part of your PATH, you can use this to add the location
# where "SRA toolkit" is installed so it can be used to download the NCBI data
# E.g.
# PATH=$PATH:/scratch1/mponce/RESEARCH/Tetrahymena_Ryerson/TOOLS/NCBI_SRA_toolkit/sratoolkit.2.9.6-1-ubuntu64/bin
# --OR--
# in Niagara/SciNet load the SRA-TOOLKIT module from CCEnv
module load CCEnv
module load StdEnv/2023 gcc/12.3
module load sra-toolkit/3.0.9


#######################################################
# check that the script is not being sourced!!!
if [[ $0 == '-bash' ]]; then
        echo "Please do not 'source' this script; ie. run it as PATHtoRACS/core/SCRIPTname arguments"
        return 1
fi

# setting preamble, detecting scripts location
scriptsDIR=`dirname $0`


# load auxiliary fns for integrity checks and message/error handling
if [[ -f $scriptsDIR/../core/auxs/auxFns.sh ]]; then
        . $scriptsDIR/../core/auxs/auxFns.sh --source-only;
else
        echo "Error auxiliary file: " $scriptsDIR/auxs/auxFns.sh "NOT found!"
        return 11
fi
#######################################################

checkTools fastq-dump gzip

# Oxytricha files:
srxFiles="SRX483016 SRX483017"

# download fastq files
echo "starting to download files using fastq-dump..."
for i in $srxFiles; do
	echo "accessing $i file ... ";
	fastq-dump -I --split-files $i
done

# generate fastq.gz files
echo "Compressing files..."
for i in *fastq; do
	echo $i ;
	gzip -c $i > $i.gz ;
done


