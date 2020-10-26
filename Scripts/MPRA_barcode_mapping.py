#!/usr/bin/env python3

"""
This script selects barcodes that will be used for downstream analysis
in RNA-seq. The perfect length requirement is not used and the Levenshtein
distance is set to the 1st percentile, not (arbitrarily) at 5
"""

import itertools
from collections import Counter, defaultdict
import Levenshtein
import numpy
import random
import argparse
import subprocess



def reverse_complement(seq):
	"""
	Return the reverse complement of a nucleotide string
	"""
	complement = {'A': 'T', 'T':'A', 'C':'G', 'G':'C'}
	
	
	rc = ''.join([complement[nt] for nt in seq[::-1]])
	return rc

def fastq_reader(filename):
	"""
	Extract sequences from FASTQ file. Each sequence takes up four lines,
	with the sequence being the second line.
	"""

	with open(filename) as infile:
		# grab lines starting at second line (1, 0-indexing), go until end of
		# file (stop = None), skip 4 lines at a time 
		line_slice = itertools.islice(infile, 1, None, 4)
		
		# convert iterator to list, strip new lines
		sequences = []
		for line in line_slice:
			sequences.append(line.strip())

	return sequences

def get_wc(filename):
	p = subprocess.Popen(['wc', '-l', filename], stdout = subprocess.PIPE,
												 stderr = subprocess.PIPE)
	out, err = p.communicate()
	lines = float(out.split()[0])
	return lines

def library_reader(filename, primer_length, rev_complement=True, format = 'csv'):
	"""
    Read in .csv or tab file of library sequences. First column is sequence name,
    second column is sequence. Trim primer sequences.
	"""

	lib = {}

	with open(filename) as infile:
		# read through first line, i.e. skips first line. Need column names in tsv
		infile.readline()
		for line in infile:
			if format == 'csv':
				name, seq = line.strip().split(',')[:2]
			elif format == 'tab':
				name, seq = line.strip().split()[:2]
			
			seq = seq.upper()
			seq = seq[primer_length:-primer_length]

			if rev_complement:
				rc_seq = reverse_complement(seq)
				lib[rc_seq] = name

			#lib[seq] = name #### could this be the bad line???? commented out

	return lib

def find_perfect_reads(reads, lib, var_length):
	"""
	Extract reads that perfectly match the library sequence
	"""
	# only take reads that perfectly match a reference sequence, get rid of length requirement
	# perfect_reads = [read for read in reads if read[:var_length] in lib]
	perfect_reads = [read for read in reads if read[-var_length:] in lib]	
	return perfect_reads

def mapping(barcodes, perfect_reads, reads, bc_loc, bc_length, var_length):

	variant_map = defaultdict(list)
	barcode_map = defaultdict(list)

	# for each barcode that passes filters, look at all reads and see what
	# it maps to
	for read in reads:
		if bc_loc == 'start':
			barcode = read[:bc_length]
		elif bc_loc == 'end':
			barcode = read[-bc_length:]

		# if barcode in filtered set
		if barcodes.get(barcode, 0) > 0:
			barcode_map[barcode].append(read)

	# in the perfect reads, keep track of how many barcodes go to each variant
	for read in perfect_reads:
		if bc_loc == 'start':
			barcode = read[:bc_length]
		elif bc_loc == 'end':
			barcode = read[-bc_length:]

		# variant = read[:var_length]
		variant = read[-var_length:]
		variant_map[variant].append(barcode)


	return [variant_map, barcode_map]


def bootstrap_levenshtein(lib, n):
	"""
	This function calculates a reference Levenshtein distribution. It randomly
	picks two sequences from the reference sequences and calculates the distance
	to get a measure of how similar the library is.
	"""

	distances = []
	# bootstrap n times
	for i in range(0, n):
		# randomly grab two sequences with replacement
		string1 = random.choice(list(lib.keys()))
		string2 = random.choice(list(lib.keys()))

		distances.append(Levenshtein.distance(string1, string2))
	
	# take cutoff at 1% percentile
	cutoff = numpy.percentile(distances, 1)

	# If the distribution consists of mainly large distances, the 1% percentile
	# will be large, so readjust the cutoff lower in this case
		 
	return cutoff

def filter_barcodes(barcode_map, cutoff, var_length, name='output.txt', lib=None):
	'''
	For each barcode, calculate the Levenshtein distance between its reads
	and if it is below cutoff (aka barcode maps to similar reads, no cross-talk)
	then keep this barcode
	'''

	final_barcodes = []
	covered_sequences = set()
	lib = set(lib.keys()) # for faster lookup

	all_dist = []

	outfile = open(name, 'w')
	headers = ['barcode', 'num_unique', 'num_reads', 'num_reads_most_common', 'most_common']
	outfile.write('\t'.join(headers)+'\n')

	for barcode in barcode_map:
		reads = barcode_map[barcode]
		# # trim off barcode (20), RE site (8), primers (15) and last 15
		# trimmed = [read[43:-15] for read in reads]
		
		# trimmed = [read[:var_length] for read in reads]
		trimmed = [read[-var_length:] for read in reads]

        # grab most common read as reference
		most_common = Counter(trimmed).most_common(1)[0][0]
		distances = [Levenshtein.distance(most_common, read) for read in set(trimmed)]
		all_dist.append(max(distances))
		# if the other reads this barcode maps to have a Levenshtein distance less than the cutoff
		if max(distances) < cutoff:
			num_unique = len(set(trimmed))
			num_reads = len(trimmed)
			num_reads_most_common = Counter(trimmed).most_common(1)[0][1]
			most_common = Counter([read[-var_length:] for read in reads]).most_common(1)[0][0]
			if most_common in lib:
				# only accept barcode if the most common read is in the library
				final_barcodes.append(barcode)
				is_reference = 1
				covered_sequences.add(most_common)
				info = [barcode, num_unique, num_reads, num_reads_most_common, most_common]
				info = map(str, info)
				outfile.write('\t'.join(info)+'\n')

	outfile.close()

	print("Percent of library represented by final barcodes:", len(covered_sequences)/ (len(lib)/2.0))

	with open('lev_dist.txt', 'w') as outfile:
		for dist in all_dist:
			outfile.write(str(dist)+'\n')

	return final_barcodes

def write_variant_results(variant_map, name, final_barcodes, lib):
	outfile = open(name, 'w')
	fields = ['variant', 'name', 'num_barcodes', 'num_unique', 'barcodes']
	outfile.write('\t'.join(fields)+'\n')

	final_barcodes = set(final_barcodes)

	for variant in variant_map:

		barcodes = variant_map[variant]
		# only keep those that are in final barcodes
		barcodes = [barcode for barcode in barcodes if barcode in final_barcodes]
		num_barcodes = len(barcodes)
		uniq_bcs = set(barcodes)
		num_unique = len(uniq_bcs)
		ref_name = lib[variant]
		info = [variant, ref_name, num_barcodes, num_unique]
		info = list(map(str, info))
		info.append(','.join(uniq_bcs))
		outfile.write('\t'.join(info)+'\n')

	outfile.close()


def check_args(args):
	if args.lib_type != 'csv' and args.lib_type != 'tab':
		raise ValueError('Please provide format type of library, either csv or tab') 

	if args.bc_loc != 'start' and args.bc_loc != 'end':
		raise ValueError('Please specify location of barcode, either start or end')


if __name__ == '__main__':

	parser = argparse.ArgumentParser('Map barcodes to sequence')
	parser.add_argument('reads_file', help='.txt file of sequences only or FASTQ file')
	parser.add_argument('lib_file', help='.csv or tab file of library')
	parser.add_argument('lib_type', help='please specify either csv or tab')
	parser.add_argument('var_length', help='length of variants', type=int)
	parser.add_argument('primer_length', help='length of primers at start and end of sequences in library file',
						type = int)
	parser.add_argument('bc_loc', help='barcode location, specify start or end')
	parser.add_argument('bc_length', help='length of barcode', type=int)
	# parser.add_argument('bc_cutoff', help='Throw out barcodes that appear less than n times', type=int)
	parser.add_argument('output_file', help='Name of output file')
	args = parser.parse_args()

	check_args(args)

	reads_file = args.reads_file
	var_length = args.var_length
	bc_loc = args.bc_loc
	bc_length = args.bc_length
	#read_length = args.read_length

	if reads_file.endswith('fastq'):
		reads = fastq_reader(reads_file)
	elif reads_file.endswith('txt'):
		reads = [line.strip() for line in open(reads_file)]
	else:
		raise Exception("Please give reads file in either FASTQ or raw one read per line txt format")

	print("Number of reads:", len(reads))

	print("Reading in library reference...")
	lib = library_reader(filename = args.lib_file, primer_length = args.primer_length,
					     rev_complement=True, format = args.lib_type)

	print("Extracting perfect reads...")
	perfect_reads = find_perfect_reads(reads = reads, lib = lib, 
									   var_length = var_length)
	print("Percent perfect:", len(perfect_reads) / float(len(reads)))

	# grab barcodes that map to a perfect sequence
	if bc_loc == 'start':
		barcodes = [read[:bc_length] for read in perfect_reads]
	elif bc_loc == 'end':
		barcodes = [read[-bc_length:] for read in perfect_reads]

	print("Number of unique barcodes for perfect reads: ", len(set(barcodes)))
	print("Filter by barcode frequency...")
	
	# Count barcodes 
	barcode_counts = dict(Counter(barcodes))

	# Throw out barcodes that appear 1 or 2 times, sequencing errors
	barcodes_clean = {x : barcode_counts[x] for x in barcode_counts if barcode_counts[x] > 2} ### changed back to greater than 2, to resolve ambiguous barcodes
	print("Number of barcodes > 2:", len(barcodes_clean))

	# barcode_cutoff = args.bc_cutoff
	# above_cutoff = {x : barcodes_clean[x] for x in barcodes_clean if barcodes_clean[x] >= barcode_cutoff}
	
	# print "Number of barcodes above cutoff:", len(above_cutoff)

	print("Mapping...")

	variant_map, barcode_map = mapping(barcodes_clean, perfect_reads, reads, bc_loc, bc_length, var_length)

	# bootstrap reference sequences to get a reference Levenshtein distribution 
	# to determine cutoff
	print("Bootstrapping reference sequences to obtain cutoff...")
	cutoff = bootstrap_levenshtein(lib, 10000)
	print("cutoff is Levenshtein distance ", cutoff)

	print("Filtering and writing results...")
	final_barcodes = filter_barcodes(barcode_map, cutoff,
									 name='barcode_statistics.txt', lib=lib, var_length = var_length)

	print("Number of final barcodes: ", len(final_barcodes))
	# write final barcodes to file
	outfile = open(args.output_file, 'w')
	for barcode in final_barcodes:
		outfile.write(barcode+'\n')
	outfile.close()

	# write variant results
	write_variant_results(variant_map, 'variant_statistics.txt', final_barcodes, lib)

