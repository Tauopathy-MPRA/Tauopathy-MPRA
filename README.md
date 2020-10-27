# Tauopathy-MPRA
Contains python mapping scripts and pre-processed data for MPRA study of Tauopathy-associated common genetic variation

## Scripts

### MPRA_barcode_mapping.py

This python3 script was used during the intermediate MPRA library construction step to identify and map DNA barcodes (appended by PCR) to unique library oligos/elements. The user must provide a previously merged paired-end sequencing .FASTQ file, or .txt file with the raw reads extracted.

Usage
```
python3 MPRA_barcode_mapping.py [FILE.*] [variant_list.*] [file type] [oligo length] [adaptor length] [bc location] [bc length] outfile.txt

mandatory arguments, in-order:

 File.*              A .FASTQ or .txt file with the raw, merged reads.
 variant_list.*      A .csv or tab file of library sequences. First column is sequence name, second column is sequence. Trim shared primer sequences
 file type           "tab" or "csv"
 olgio length        numeric length of library element/oligo (without shared adaptor sequences)
 adaptor length      numeric length of shared adaptor sequence
 bc location         location in read of barcode in relation to the oligo. "start" or "end"
 bc length           numeric length of barcodes
 outfile.txt         name of output file

Example: #used in this analysis
python3 MPRA_barcode_mapping.py BC_map_AD_PSP_2.fastq AD_PSP_variants_list_2.tsv tab 162 24 start 20 barcode_statistics.txt

```

### MPRA_BC_counter.py
This python3 script counts the number of perfect barcode reads mapped to each unique MPRA library oligo. User runs the script in subdirectory containing demultiplexed .FASTQ files. User must also provide a variant - barcode lookup table. The merge_counts.txt output file contains barcode count data that is used for subsequent downstream analysis. 

Usage

```
./MPRA_BC_counter.py [variant_statistics.txt]

mandatory arguments:
 variant_statistics.txt    tab-separated file containing barcode-variant lookup table. must have column headings "name" and "barcodes". "Name" is the oligo name, and "barcodes" corresponds to a comma separated list of each unique barcode
 
 Output
 merged_counts.txt         long-format tab-separated file. Columns correspond with unique barcode, corresponding variant name, and counts from each unique .FASTQ/replicate in the subdirectory
 count_assignments.txt      tab-separated file showing counts of the assigned (mapped) and unassigned reads for each .FASTQ/replicate.
 
 ```

## Raw Data
Raw data from this study are provided in the gzipped tar file: raw_data.tar.gz. This contains two subdirectories corresponding to data from MPRA 1 and 2 study stages. Here, raw data is preprocessed, and refers to mapped barcode counts (merged_counts.txt file generated using the MPRA_BC_counter.py script), useful for downstream analysis (see above for file format). Variant metadata files are also provided.

Meta data columns:
 group          allele: "ref" or "alt"
 ID             unique MPRA ID# for a given variant
 rsID           variant rsID
 chr            chromosome
 pos            position (hg19)
 A0             hg19 reference allele
 A1             hg19 alternate allele
 RegulomeDB     RegulomeDB score 
 variant        162 bp unique oligo
 name           unique name for each oligo/element. Format: group_ID_type
 
 -note on types: For MPRA 2, under name there are additional variant type annotations: 
    "new" refers to a newly tested variant. 
    "missing" refers to a variant re-tested in MPRA 2 that had dropped out of MPRA 1. 
    "negative" refers to 50 non-significant MPRA 1 variants retested in MPRA 2. 
    "fdr" refers to 212 exactly replicated MPRA 1 SigVars. 
    "RC" refers to MPRA 1 SigVars in reverse complement orientation. 
    "upstream" are SigVars with the allele placed in the lower third of genomic context rather than the middle (standard oligo design). 
    "reverse" are SigVars in reverse orientation.
    
-- Raw, unmapped data (unprocessed .FASTQ files) for this study are uploaded to SRA, accession: 

-- If using scripts or raw data from this study, please cite:
