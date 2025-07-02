#!/usr/bin/env python3

# RACS: InterGenic Regions -- Python main driver script
# Requires: utils_RACS_IGR.py
# Usage: python interGeneRegions.py inputFile refFILE.gff3 outFile

from utils_RACS_IGR import (
    parse_args, read_data, read_ref_table, reshape_table,
    dump_data, save_data, scfld_cap, warning_scfld
)

def main():
    print(" RACS v1.0 (2018/2019) -- InterGenic Regions Finder")
    print("----------------------------------------------------")

    input_file, ref_file, output_file = parse_args()
    input_data = read_data(input_file)
    ref_table = read_ref_table(ref_file)
    sorted_data = reshape_table(input_data)

    nbr_entries = len(sorted_data)
    lstscfld, lstregion1, lstregion2, lstSize = [], [], [], []
    SCFLD = ""
    endRegion = None

    for i in range(nbr_entries):
        scaffold = sorted_data.loc[i, 'lscaffold']
        begReg = sorted_data.loc[i, 'lbregion']
        endReg = sorted_data.loc[i, 'leregion']
        print(i + 1, scaffold, begReg, endReg)

        if SCFLD == scaffold:
            if i > 0 and begReg > sorted_data.loc[i - 1, 'leregion']:
                beginRegion = begReg
                myScfld = dump_data(SCFLD, endRegion, beginRegion, 1)
                lstscfld.append(myScfld[0])
                lstSize.append(myScfld[3])
                lstregion1.append(endRegion + 1)
                lstregion2.append(beginRegion - 1)
                endRegion = endReg
            else:
                warning_scfld(scaffold, begReg, endReg, lstscfld, lstregion1, lstregion2)
                endRegion = endReg  # Still update endRegion
        else:
            if SCFLD != "":
                beginRegion = begReg
                scfCAP = scfld_cap(SCFLD, ref_table)
                if scfCAP:
                    if endRegion < scfCAP:
                        myScfld = dump_data(SCFLD, endRegion, scfCAP + 1, 1)
                        lstscfld.append(myScfld[0])
                        lstSize.append(myScfld[3])
                        lstregion1.append(endRegion + 1)
                        lstregion2.append(scfCAP)
                    else:
                        warning_scfld(scaffold, begReg, endReg, lstscfld, lstregion1, lstregion2)

                    myScfld = dump_data(scaffold, 0, beginRegion, 1)
                    lstscfld.append(myScfld[0])
                    lstSize.append(myScfld[3])
                    lstregion1.append(1)
                    lstregion2.append(beginRegion - 1)
                    endRegion = endReg
                else:
                    print("Potential data inconsistency issue detected... please verify your data integrity...\n")
            else:
                # first scaffold
                beginRegion = begReg
                endRegion = endReg
                myScfld = dump_data(scaffold, 0, beginRegion, 1)
                lstscfld.append(myScfld[0])
                lstSize.append(myScfld[3])
                lstregion1.append(1)
                lstregion2.append(beginRegion - 1)

            SCFLD = scaffold

        if len(lstregion1) != len(lstregion2):
            raise RuntimeError("Mismatch in region list lengths")

    # Handle last case
    scfCAP = scfld_cap(scaffold, ref_table)
    if scfCAP:
        myScfld = dump_data(scaffold, endRegion, scfCAP + 1, 1)
        lstscfld.append(myScfld[0])
        lstSize.append(myScfld[3])
        lstregion1.append(endRegion + 1)
        lstregion2.append(scfCAP)

    # Output final intergenic regions
    import pandas as pd
    interGenes = pd.DataFrame({
        "scaffold": lstscfld,
        "beggining": lstregion1,
        "end": lstregion2,
        "size": lstSize
    })

    save_data(interGenes, output_file)

if __name__ == "__main__":
    main()