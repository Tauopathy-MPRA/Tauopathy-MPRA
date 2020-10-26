# Tauopathy-MPRA
Contains python mapping scripts and pre-processed data for MPRA study of Tauopathy-associated common genetic variation

## Scripts

### MPRA_barcode_mapping.py

This python3 script was used during the intermediate MPRA library construction step to identify and map DNA barcodes (appended by PCR) to unique library oligos/elements. The user must provide a previously merged paired-end sequencing .FASTQ file, or .txt file with the raw reads extracted.

Usage
```
python3 MPRA_barcode_mapping.py [FILE.*] [variant_list.*] [file type] [oligo length] [adaptor length] [bc location] [bc length] outfile.txt

mandatory arguments, in-order:

 File.*              A .FASTQ or .txt file with the raw, merged sequence reads.
 variant_list.*      A .csv or tab file of library sequences. First column is sequence name, second column is sequence. Trim primer sequences
 file type           "tab" or "csv"
 olgio length        numeric length of library element/oligo (without shared adaptor sequence)
 adaptor length      numeric length of shared primer sequence on flanking oligo
 bc location         location in read of barcode in relation to the oligo. "start" or "end"
 bc length           numeric, length of unique barcode
 outfile.txt         name of output file

Example: #used in this analysis
python3 MPRA_barcode_mapping.py BC_map_AD_PSP_2.fastq AD_PSP_variants_list_2.tsv tab 162 24 start 20 barcode_statistics.txt

```

### MPRA_BC_counter.py
This python3 script counts the number of perfect barcodes mapped to each unique MPRA library oligo. User runs script in directory containing demultiplexed .FASTQ files. User must also provide a variant - barcode lookup table. The merge_counts.txt output file contains barcode count data that is used for subsequent downstream analysis. 

Usage

```
./MPRA_BC_counter.py [variant_statistics.txt]

mandatory arguments:
 variant_statistics.txt    tab-separated file containing barcode-variant lookup table. must have column headings "name" and "barcodes". "Name" is the oligo name, and "barcodes" corresponds to a comma separated list of each unique barcode
 
 Output
 merged_counts.txt         long-format tab-separated file. Columns correspond with unique barcode, its corresponding variant, and counts from each unique .FASTQ/replicate in the subdirectory
 count_assignments.txt      tab-separated file showing counts of the assigned (mapped) and unassigned reads for each .FASTQ/replicate.
 
 ```



