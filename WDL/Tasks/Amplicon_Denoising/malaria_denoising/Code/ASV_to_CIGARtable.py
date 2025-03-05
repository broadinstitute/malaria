#!/usr/bin/env python

import os
import sys
import argparse
#import subprocess
#import csv
#import gzip
#import pandas as pd

#from Bio import SeqIO

import asv_to_cigar as ac

def main():
	### LOAD ARGUMENTS

	#Parse command line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("--reference_amplicons", type=str, help="Path to reference amplicons")
	parser.add_argument("--seqtab", type=str, help="Path to seqtab file")
	parser.add_argument("--ASVSeqs", type=str, help="Path to fasta file")
	parser.add_argument("--ASVTable", type=str, help="Path to table file")
	parser.add_argument("--min_reads", type=int, help="Minimum reads threshold")
	parser.add_argument("--min_samples", type=int, help="Minimum samples threshold")
	parser.add_argument("--max_snv_dist", type=int, help="Maximum SNV distance")
	parser.add_argument("--max_indel_dist", type=int, help="Maximum indel distance")
	parser.add_argument("--include_failed", action="store_true", default=False, help="Include failed samples in analysis")
	parser.add_argument("--exclude_bimeras", action="store_true", default=False, help="Exclude bimeras from analysis")
	parser.add_argument("--polyN", type=int, help="Minimum homopolymer length")

	args = parser.parse_args()				

	### PREPARE OUTPUT DIRECTORIES
	#Generate the path to the Results directory.
	global res_dir
	res_dir = os.path.abspath(os.path.join("Results"))

	#ASV to CIGAR
	#Convert ASVs from DADA2 pipeline to pseudo-CIGAR strings.
	print("Converting ASVs to CIGARs")

	#Load arguments
	path_to_amp_db = args.reference_amplicons
	path_to_seqtab = args.seqtab
	path_to_fasta = args.ASVSeqs
	path_to_table = args.ASVTable
	min_reads = args.min_reads
	min_samples = args.min_samples
	max_snv_dist = args.max_snv_dist
	max_indel_dist = args.max_indel_dist
	include_failed = args.include_failed
	exclude_bimeras = args.exclude_bimeras
	polyN = args.polyN

	#Output files
	path_to_alignments = os.path.join(res_dir, "Alignments") #Directory to store ASV alignment files
	path_to_out = os.path.join(res_dir, "CIGARVariants_Bfilter.out.tsv") #Output seqtab tsv file with amplicon/variant counts
	path_asv_to_cigar = os.path.join(res_dir, "ASV_to_CIGAR.out.txt") #Output file for ASV -> CIGAR string table 
	path_to_zero_read_samples = os.path.join(res_dir, "ZeroReadsSampleList.txt") #Output file for 

	print(f"INFO: Loading {path_to_amp_db}")
	amplicons = ac.parse_amp_db(path_to_amp_db)
	if not amplicons:
		print(f"ERROR: No amplicons in {path_to_amp_db}")
		sys.exit(1)

	#TO DELETE? THIS MASK CURRENTLY NOT BEING USED.
	#if os.path.exists("amp_mask.txt"):
	#	print(f"INFO: Loading amp_mask.txt")
	#	mask = ac.parse_dustmasker("amp_mask.txt")
	#else:
	#	print(f"INFO: No mask data specified.")
	#	mask = {}
	mask = {} #Provisonal for ac.parse_alignments, remove if the mask is command is deleted. 

	print(f"INFO: Loading {path_to_fasta}")
	asvs = ac.get_asv_seqs(path_to_fasta)
	if not asvs:
		print(f"ERROR: No ASV sequences in {path_to_fasta}")
		sys.exit(1)

	print(f"INFO: Parsing {path_to_table} with total reads >= {min_reads}, samples >= {min_samples}, snv_dist <= {max_snv_dist}, indel_dist <= {max_indel_dist}")

	if include_failed:
		print("WARNING: Including ASVs that failed post-DADA2 filters! This is not recommended.")
	else:
		print("INFO: Excluding ASVs that failed post-DADA2 filters.")

	if exclude_bimeras:
		print("INFO: Excluding ASVs that DADA2 marked as bimeras.")

	bins = ac.parse_asv_table(path_to_table, min_reads=min_reads, min_samples=min_samples, max_snv_dist=max_snv_dist, max_indel_dist=max_indel_dist, include_failed=include_failed, exclude_bimeras=exclude_bimeras) #This function only matches to the first strain.
	if not bins:
		print(f"ERROR: No useable data in {path_to_table}")
		sys.exit(1)

	print(f"INFO: Writing amplicon fasta files to {path_to_alignments}")
	ac.write_amplicon_fastas(asvs, bins, amplicons, outdir=path_to_alignments)

	print("INFO: Running MUSCLE aligner on amplicon fasta files. Please wait...")
	ac.run_muscle(bins, outdir=path_to_alignments)

	print("INFO: Parsing alignments to CIGAR strings")
	cigars = ac.parse_alignments(bins, mask=mask, min_homopolymer_length=polyN, outdir=path_to_alignments)
	if not cigars:
		print("ERROR: could not determine CIGAR strings")
		sys.exit(1)

	if path_asv_to_cigar:
		ac.write_cigar_strings(cigars, path_asv_to_cigar)
		print(f"INFO: Wrote ASV->CIGAR table to {path_asv_to_cigar}")

	print(f"INFO: Converting DADA2 seqtab file {path_to_seqtab} to {path_to_out}")
	if ac.convert_seqtab(path_to_seqtab, cigars, path_to_out):
		print("INFO: Completed conversion of seqtab to CIGAR variants successfully!")

		if ac.get_zero_reads_samples(path_to_out, path_to_zero_read_samples):
			print("INFO: Obtained samples with zero reads successfully!")
		
if __name__ == "__main__":
	main()
