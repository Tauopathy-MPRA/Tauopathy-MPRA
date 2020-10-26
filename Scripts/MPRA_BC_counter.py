#!/usr/bin/env python3

import pandas as pd
import os
import sys
from itertools import islice
from functools import reduce

def make_barcode_dict(variant_list):
    df = pd.read_csv(variant_list, sep='\t', usecols=["name", "barcodes"])
    barcode_dict = {}
    barcode_dup = {}
    for i in df.itertuples():
        if type(i[2]) == float:
            continue
        else:
            for j in i[2].split(","):
                if j in barcode_dict:
                    barcode_dup[j] = i[1]
                else:
                    barcode_dict[j] = (i[1], 0)
    entriesToRemove = barcode_dup.keys()
    for k in entriesToRemove:
        barcode_dict.pop(k, None)
    return barcode_dict

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])

def read_mapping(reads_rc, barcode_dict, inputdir, name):
    assigned_count = 0
    unassigned_count = 0
    #unassigned_reads = []
    update_dict = barcode_dict.copy()  # shallow copy barcode_dict to a dictionary for editing
    for i in reads_rc:
        if i in update_dict.keys():
            update_dict[i] = (update_dict[i][0], update_dict[i][1] + 1)
            assigned_count += 1
        else:
            #unassigned_reads.append(i)
            unassigned_count += 1
    #with open(inputdir + name + "_unassigned_reads", "w") as g:  # writes the unassigned reads to a text file
        #g.writelines("%s\n" % line for line in unassigned_reads)
        #g.close()
    return update_dict, assigned_count, unassigned_count

def read_process(filelist, inputdir, barcode_dict, num_files):
    files_processed = 0
    master_count_dict = {}
    master_assigned_count = {}
    master_unassigned_count = {}
    for i in filelist:  # iterates through all files in the working directory
        if i.endswith(".fastq"): # only looks at fastq files
            name = i.split("_")[0]
            reads_rc = []
            with open(inputdir+i, "r") as f:
                for line in islice(f, 1, None, 4):  # gets every fourth line to extract only reads
                    linestrip = line.strip()  # removes whitespace
                    readstrip = linestrip[:20]  # gets first 20bp of 26bp reads
                    reads_rc.append(reverse_complement(readstrip))  # gets reverse complement of the reads, to match to barcode directory
            master_count_dict[name], master_assigned_count[name], master_unassigned_count[name] = read_mapping(reads_rc, barcode_dict, inputdir, name)
            f.close()
            files_processed += 1
            print(files_processed / num_files * 100, '% complete ')
    return master_count_dict, master_assigned_count, master_unassigned_count

def master_count_df(master_count_dict):
    master_df_list = []
    for i in master_count_dict:
        df = pd.DataFrame.from_dict(master_count_dict[i], orient='index', columns=['variant', i])
        df.index.name = 'barcode'
        df.reset_index(inplace=True)
        master_df_list.append(df)
    df_merged = reduce(lambda left, right: pd.merge(left, right, on=['barcode', 'variant'], how='inner', suffixes=(None, None)), master_df_list) 
    pd.DataFrame.to_csv(df_merged, inputdir + "merged_counts.txt", sep= '\t', index = False)

def master_assignment_df(master_assigned_count, master_unassigned_count):
    df1 = pd.DataFrame.from_records(master_assigned_count, index=['assigned_counts'])
    df2 = pd.DataFrame.from_records(master_unassigned_count, index=['unassigned_counts'])
    df3 = pd.concat(objs=[df1,df2], axis=0)
    df3.index.name = 'count_type'
    df3.reset_index(inplace=True)
    pd.DataFrame.to_csv(df3, inputdir + "count_assignments.txt", sep='\t', index= False)


if __name__ == "__main__":

    print("reading directories")
    inputdir = sys.path[0] + "/"
    variant_list = inputdir + sys.argv[1]
    filelist = os.listdir(inputdir)
    num_files = 0
    for i in filelist:
        if i.endswith(".fastq"):
            num_files += 1

    print("processing variant list")
    barcode_dict = make_barcode_dict(variant_list)
    print("processing FastQ files")
    master_count_dict, master_assigned_count, master_unassigned_count = read_process(filelist, inputdir, barcode_dict, num_files)
    print("printing combined counts files")
    master_count_df(master_count_dict)
    master_assignment_df(master_assigned_count, master_unassigned_count)
