#!/usr/bin/env python

import os
import sys
import argparse
import subprocess
import csv

from Bio import SeqIO

import ci_detection as ci

def main():
	"""
	Implementation of the amplicon contamination detection pipeline for use and 
	distribution in TERRA v=1.0

	Usage: python Code/CI_TerraPipeline.py 
	
	Returns:
	None. However, see the functions documentation.
	"""

	### LOAD ARGUMENTS

	#Parse command line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('--adaptors_csv', '-a', type=str, help="CSV file containing fastq1s and fastq2s")
	parser.add_argument('--barcodes', '-b', type=str, help="List of barcodes")
	parser.add_argument('--found_files', '-l', type=str, help="List of found files")
	parser.add_argument('--distance', '-d', type=int, default=0, help="Distance for CI") 

	args = parser.parse_args()

	adaptor_rem1s = []
	adaptor_rem2s = []
	with open(args.adaptors_csv, newline='') as csvfile:
		reader = csv.DictReader(csvfile)
		for row in reader:
			adaptor_rem1s.append(row['adaptor_rem1'])
			adaptor_rem2s.append(row['adaptor_rem2'])

	barcodes = args.barcodes
	found_files = args.found_files
	distance = args.distance
	print(f"Adaptor_rem1s: {adaptor_rem1s}")
	print(f"Adaptor_rem2s: {adaptor_rem2s}")
				
	### PREPARE OUTPUT DIRECTORIES
	#Generate the path to the Results directory.
	global res_dir 
	global rep_dir
	res_dir = os.path.abspath(os.path.join("Results"))
	rep_dir = os.path.abspath(os.path.join("Report"))
	os.mkdir(res_dir)
	os.mkdir(rep_dir)

	#Merge forward and reverse reads with bbmerge. Only for reads that overlap and
	#Process merge report and generate RStudio plots
	ci.flush_dir(res_dir, "Merge")
	with open(args.found_files, 'r') as f:
		samples = [line.strip() for line in f]

	for sample in samples:
		# Find matches in adaptor_rem1s and adaptor_rem2s
		match1 = next((f for f in adaptor_rem1s if sample in f), None)
		match2 = next((f for f in adaptor_rem2s if sample in f), None)

		print(f"Processing sample {sample}")
		print(f"Match1: {match1}")
		print(f"Match2: {match2}")

		if match1 and match2:
			ci.mergereads(sample, match1, match2, res_dir, "Merge")
		else:
			print(f"No matches found for sample {sample} in adaptor_rem1s and adaptor_rem2s")

	ci.flush_dir(rep_dir, "Merge")
	meta = open(os.path.join(res_dir, "Merge", "merge_meta.tsv"), 'r')
	samples = meta.readlines()
	
	print("Extracting fields from bbmerge reports")
	for sample in samples:
		slist = sample.split()
		ci.extract_bbmergefields(slist[0], slist[1], slist[3], barcodes, res_dir, rep_dir, "Merge", distance)
		
if __name__ == "__main__":
	main()
