# Tauopathy-MPRA
Contains python mapping scripts and preprocessed data for MPRA study of Tauopathy-associated common genetic variation

## Scripts
### MPRA_barcode_mapping.py

This python3 script was used during the intermediate MPRA library construction step to identify and map DNA barcodes appended by PCR to unique library oligos/elements. The user must provide a previously merged paired-end sequencing .FASTQ file, or .txt file with the raw reads extracted.

Example Usage
```
python3 MPRA_barcode_mapping.py [FILE.*/] [variant_list.*] [file type] [oligo length] [adaptor length] [bc location] [bc length] outfile.txt

mandatory arguments, in-order:
File.*  A .FASTQ or .TXT file with the raw, merged sequence reads.
variant_list.*  A .csv or tab file of library sequences. First column is sequence name, second column is sequence. Trim primer sequences
file type "tab" or "csv"
olgio length  numeric length of library element/oligo (without shared adaptor sequence)
adaptor length  numeric length of shared primer sequence on flanking oligo
bc location location in read of barcode in relation to the oligo. "start" or "end"
bc length numeric, length of unique barcode
outfile.txt name of output file

Example: #used in this analysis
python3 MPRA_barcode_mapping.py BC_map_AD_PSP_2.fastq AD_PSP_variants_list_2.tsv tab 162 24 start 20 barcode_statistics.txt

```


python3 barcode_mapping_less_strict.py BC_map_AD_PSP_2.fastq AD_PSP_variants_list_2.tsv tab 162 24 start 20 barcode_statistics.txt
