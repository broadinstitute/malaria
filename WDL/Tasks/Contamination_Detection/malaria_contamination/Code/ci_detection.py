"""
Combinatonrial indices detection scripts
"""
import os, fnmatch
import subprocess
import sys
import shutil
import glob

from Bio import SeqIO
from Bio.Seq import Seq

#GENERAL USE SCRIPTS SECTION
def flush_dir(parent_dir, *dirnames):
	"""
	Remove all files and subdirectories from a directory, and create a new empty
	directory with the same name. Multiple subdirectiories may be provided.

	Args:
	parent_dir (str): The path of the parent directory.
	dirname (str): The name of the directory to be flushed.

	Returns:
	None
	"""

	dirpath = os.path.join(parent_dir, *dirnames)
	shutil.rmtree(dirpath, ignore_errors=True)
	os.makedirs(dirpath)
	return ()

def create_meta(path_to_fq, parent_dir, dirname, filename, pattern_fw, pattern_rv):
	"""
	Creates a metadata file with the list of files to process and their paths.

	Args:
	path_to_fq: string, path to the directory containing the input fastq files
	parent_dir: string, path to the parent directory where the output files will be stored
	dirname: string, name of the subdirectory where the output files will be stored
	filename: string, name of the output metadata file
	pattern_fw: string, pattern to match forward reads
	pattern_rv: string, pattern to match reverse reads

	Returns:
	None

	Example usage:

	create_meta('/path/to/input/fastq/', '/path/to/output/', 'subdir', 'metadata.tsv', '*_R1.fastq', '*_R2.fastq')

	This will search for all files in /path/to/input/fastq/ that end with '_R1.fastq' and '_R2.fastq',
	create a metadata file named 'metadata.tsv', and store it in a subdirectory named 'subdir' within 
	/path/to/output/.
	"""

	filelist = os.listdir(path_to_fq)

	dirpath = os.path.join(parent_dir, dirname) 
	outfile = os.path.join(dirpath, filename)

	meta_df = []
	for file_fw in glob.glob(os.path.join(path_to_fq, pattern_fw)):
		sampleid = os.path.basename(file_fw).split(pattern_fw[1:])[0]
		file_rv = os.path.join(path_to_fq, sampleid + pattern_rv[1:])
		if os.path.isfile(file_rv):
			meta_df.append((sampleid, file_fw, file_rv))

	with open(outfile, 'w') as f:
		for row in meta_df:
			f.write('\t'.join(row) + '\n')

	print(f"Meta file generated at location {outfile}")
	return()

def mergereads(sampleid, fileF, fileR, res_dir, subdir):
	"""
	This function uses bbmerge.sh to merge paired-end reads from two fastq files
	(fileF and fileR) into a single fastq file. It also generates two other fastq
	files of unmerged reads. The output files are saved in the specified res_dir
	and subdir directory paths. The function also creates a metadata file 
	(merge_meta.tsv) containing the sample ID, output file name, and standard 
	output and error logs. If either the forward or reverse fastq files are not 
	found, the function exits with an error message. This functions is optimized 
	for reads that are at the least 200bp long and amplicons 100 bp long or longer.
	Merging shorter reads will require chaning this parameters in the config file.
	
	Args:
	sampleid: a string representing the sample identifier.
	fileF: a string representing the file path of the forward reads.
	fileR: a string representing the file path of the reverse reads.
	res_dir: a string representing the directory path of the results.
	subdir: a string representing the subdirectory path of the results.

	Returns: None

	Example usage:

	mergereads("sample1", "/path/to/forward.fastq", "/path/to/reverse.fastq", "/path/to/results", "subdirectory")
	"""

	if os.path.isfile(fileF) and os.path.isfile(fileR):	
		file_nameout = os.path.join(res_dir, subdir, f"{sampleid}_stdout.txt")
		file_nameerr = os.path.join(res_dir, subdir, f"{sampleid}_stderr.txt")
		output_file_path = os.path.join(res_dir, subdir, f"{sampleid}_merged.fastq")
		output_unmerged_f_path = os.path.join(res_dir, subdir, f"{sampleid}_unmergedf.fastq")
		output_unmerged_r_path = os.path.join(res_dir, subdir, f"{sampleid}_unmergedr.fastq")
		meta_file_path = os.path.join(res_dir, subdir, f"merge_meta.tsv")

		original_stdout = sys.stdout
		original_stderr = sys.stderr

		sys.stdout = open(file_nameout, "w")
		sys.stderr = open(file_nameerr, "w")

		with open(meta_file_path, "a") as meta_file:
			meta_file.write(f"{sampleid}\t{output_file_path}\t{file_nameout}\t{file_nameerr}\n")

# OPTIONAL FOR  CONCACTENATING READS
#		cmd = ['fuse.sh', 
#			f'in1={fileF}',
#			f'in2={fileR}',
#			f'out={output_file_path}',
#			f'fusepairs=t',
#			f'pad=10']
#

		cmd = ['bbmerge.sh', 
			f'in1={fileF}',
			f'in2={fileR}',
			f'out={output_file_path}']

		print(cmd)
		proc = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
		proc.wait()

		sys.stdout = original_stdout
		sys.stderr = original_stderr
	else:
		sys.exit('BBmerge halted : one or both of the fastq files not found! Exiting..')
	return()

def extract_bbmergefields(sampleid, mergefile, bbreportfile, path_to_barcodes, res_dir, rep_dir, subdir, distance):
	"""
	Extracts relevant data from a bbmerge report file and saves it to a tab-separated file.

	Args:
	sampleid: the ID of the sample being processed
	mergefile: the path to the file with the merged reads
	bbreportfile: the path to the bbmerge report file
	path_to_barcodes: the path to a csv file with the sample_id,Forward,Reverse, where Forward and Reverse are columns with the barcodes for the sample
	res_dir: the path to the main results directory
	rep_dir: the path to the reports directory within the results directory
	subdir: the name of the subdirectory within the results directory where output files should be written

	Returns:
	None

	Example Usage:

	extract_bbmergefields("Sample1", "/path/to/bbmerge.fastq", "/path/to/bbmerge_report.txt", "path/to/barcodes_match.csv", "/path/to/results", "/path/to/reports", "bbmerge")
	"""

	if os.path.isfile(bbreportfile) and os.path.isfile(mergefile):				
		bbmergedata = {}
		bbmergedata['sampleid'] = sampleid

		path_to_program = os.path.join("/", "Code/", "runBBMergecontamination.R")
		platform = '--terra'

		cmd = ['Rscript', path_to_program,
		'-p', f'{mergefile}',
		'-d', os.path.join(rep_dir, subdir),
		'-b', path_to_barcodes,
		'-t', f'{distance}']

		cmd.append(f'{platform}')

		print(cmd)
		proc = subprocess.Popen(cmd)
		proc.wait()
	else:
		sys.exit('Extract bbmerge report halted : bbmerge report file not found! Exiting..')
	return()
