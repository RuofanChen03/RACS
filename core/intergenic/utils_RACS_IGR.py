import pandas as pd
import re
import os
import sys
import argparse

def reshape_table(origDATA):
    lscaffold = []
    lbregion = []
    leregion = []

    # Reshape the table to FIX order
    for i in origDATA['geneSCFFLD']:
        scaffold, region = i.split(':')
        lscaffold.append(scaffold)
        bregion, eregion = map(int, region.split('-'))
        lbregion.append(bregion)
        leregion.append(eregion)

    tmpTable = pd.DataFrame({
        'lscaffold': lscaffold,
        'lbregion': lbregion,
        'leregion': leregion
    })

    # Sort using numeric part of the scaffold
    scfSTRG = [int(re.sub(r'\D', '', s)) if re.search(r'\d', s) else 0 for s in lscaffold]
    tmpTable['scfSTRG'] = scfSTRG
    sortedTABLE = tmpTable.sort_values(by='scfSTRG').drop(columns=['scfSTRG']).reset_index(drop=True)

    return sortedTABLE

def dump_data(scaffold, region1, region2, flag):
    r1 = int(region1) + 1
    if region2 != 'xxx':
        r2 = int(region2) - 1
        geneSz = r2 - r1 + 1
    else:
        r2 = region2
        geneSz = 'XXX'

    strng0 = f"{scaffold}:{r1}-{r2}"
    strng1 = f"{strng0}\t{geneSz}"

    if flag == 1:
        print(strng1)

    return [strng0, r1, r2, geneSz]

def save_data(data, file_name):
    # Assuming `data` is a pandas DataFrame
    data.to_csv(file_name, sep='\t', index=False, quoting=3)  # quoting=3 means csv.QUOTE_NONE

def err_msg_fn(*args):
    print("RACS: intergenic region determination script")
    print("----------------------------------------------------")
    print("Attention! This script requires 3 arguments!")
    print("\t i: input file where to read the combined tables from")
    print("\t ii: reference *gff3* genome file for the organism, e.g., 'T_thermophila_June2014.sorted.gff3'")
    print("\t iii: name of the generated intergenic files, e.g., 'interGENs.csv'")
    print("----------------------------------------------------\n")
    print(''.join(args))
    sys.exit(1)

def check_file(file_name):
    if not os.path.exists(file_name):
        err_msg_fn(f"Error: '{file_name}' NOT found!")

def parse_args(default_ref_file="DATA/T_thermophila_June2014.sorted.gff3", default_out_name="interGENs.csv"):
    parser = argparse.ArgumentParser(description="RACS intergenic region pipeline")
    parser.add_argument("inputFile", help="Input file to read combined tables from")
    parser.add_argument("refFile", nargs="?", default=None, help="Reference gff3 file")
    parser.add_argument("outFile", nargs="?", default=None, help="Output file name")

    args = parser.parse_args()
    input_file = args.inputFile
    check_file(input_file)

    if args.refFile:
        ref_file = args.refFile
        check_file(ref_file)
        out_file = args.outFile if args.outFile else default_out_name
        if not args.outFile:
            print(f"Assuming default name for output file... '{default_out_name}'")
    else:
        err_msg_fn("Error: this script requires at least 2 arguments!")

    return input_file, ref_file, out_file

def read_data(filename):
    input_data = pd.read_csv(filename, header=0, sep='')  # assuming comma delimiter like default R `read.csv`
    input_data.rename(columns={input_data.columns[0]: "geneSCFFLD"}, inplace=True)

    return input_data

def read_ref_table(ref_file, kwrd='contig'):
    ref_table_orig = pd.read_csv(ref_file, header=None, sep='\t')

    # Filter rows where column 2 (V3 in R) contains the keyword (partial match)
    matches = ref_table_orig[ref_table_orig[2].str.contains(kwrd, na=False)]
    ref_table = matches[[0, 4]].copy()  # columns V1 and V5 in R (0-based in Python)

    orig_sz = len(ref_table)
    ref_table = ref_table.drop_duplicates()
    
    print(f"Original records in ref. table: {orig_sz}")
    print(f"After eliminating possible duplicates: {len(ref_table)}")

    ref_table.columns = ["scaffold", "supercontig"]
    return ref_table

def scfld_cap(scfld, ref_table):
    result = ref_table[ref_table['scaffold'] == scfld]['supercontig']
    return result.values[0] if not result.empty else None

def warning_scfld(scfld, beg_reg, end_reg, lstscfld, lstregion1, lstregion2):
    print(">>>>>>>>>>>>>>>> SUSPICIOUS OVERLYING REGIONS!!! <<<<<<<<<<<<<<<<")
    print(scfld, beg_reg, end_reg)
    try:
        print(lstscfld[-2], lstregion1[-2], lstregion2[-2])
    except IndexError:
        print("Warning: Not enough entries in the lists to access previous elements.")