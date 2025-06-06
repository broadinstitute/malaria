import argparse
import pandas as pd
import sys
from collections import defaultdict 
import os
import subprocess 
import re

def main():
    """
    Loads at all relevant input and intermediate files and checks that all header names match. This makes 
    the pipeline more robust to crashing from user input.   

    Usage: python Code/check_input_headers.py -i [panel_info]

    Returns:
    None.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="Input path to amplicon panel_info file.", required=True)
    parser.add_argument('-rg', '--ref_genome', help="Input path for reference genome.", required=True)
    parser.add_argument('--ref_amplicons_provided', help="Boolean if reference panel was provided or auto-generated.", action='store_true')

    args = parser.parse_args()

    panel_info = args.input
    ref_genome = args.ref_genome
    ref_amplicons_provided = args.ref_amplicons_provided

    # All of these files are either user-provided OR generated by the program. Failure to have these files means there is an upstream issue.
    panel_info_bed = "amplicon_panel.bed"
    ref_file = "reference.fasta"
    fw_file = "primers_fw.fasta"
    rv_file = "primers_rv.fasta"
    markers_table_file = "markers_table.csv"

    # Check 1 - If reference_amplicons_provided, compare it to reference that is generated from panel_info. Generate markers_table accordingly.
    if ref_amplicons_provided:
        # Reference amplicons may have different names from what could be autogenerated
        with open("reference.generated.fasta", 'w') as fd:
            subprocess.run(["bedtools", "getfasta", "-fi", ref_genome, "-bed", panel_info_bed], stdout=fd)
        # Compare autogenerated fasta to reference fasta
        ref_generated_dict, ref_generated_header_to_index = createHeaderSeqDict("reference.generated.fasta")
        ref_provided_dict, ref_provided_header_to_index = createHeaderSeqDict("reference.fasta")
        
        seqs_generated = set(ref_generated_dict.values())
        seqs_provided = set(ref_provided_dict.values())

        seqs_generated_1clip = {s[1:] for s in seqs_generated}
        seqs_provided_1clip = {s[1:] for s in seqs_provided}

        seqs_to_use = seqs_generated        
        if ref_generated_dict != ref_provided_dict: 
            # Check 1.1 - The sequences must be equal or 1 bp clip away.
            if (seqs_generated == seqs_provided):
                print("The sequences between panel_info + reference_genome and reference_amplicons are equal.")
            # Special cases due to 1-nucleotide clipping of certain amplicon reference sequences
            elif (seqs_generated_1clip == seqs_provided):
                print("The generated sequences have 1 extra base pair relative to the provided sequences. Using provided sequences...")
                seqs_to_use = seqs_provided
            elif (seqs_provided_1clip == seqs_generated):
                print("The provided sequences have 1 extra base pair relative to the generated sequences. Using generated sequences...")
                seqs_to_use = seqs_generated
            # Unequal references -> Incorrect panel or sequence information
            else:
                print("The sequences between the panel_info + reference_genome and the reference_amplicons must be equal.")
                print("Please provide the correct reference_genome and sequence positions in panel_info, or the correct sequences in reference_amplicons.")
                sys.exit()
            
            # Check 1.2b - If the headers are not equal, correct
            equal_headers = set(ref_generated_dict.keys()) == set(ref_provided_dict.keys())
            if equal_headers:
                print("The headers are equal.")
            else:
                print("The headers are not equal. Replacing auto-generated FASTA headers with provided headers from reference_amplicons...")

            # Create sequence-to-header mappings
            seq_to_generated_header = {v: k for k, v in ref_generated_dict.items()}
            seq_to_provided_header = {v: k for k, v in ref_provided_dict.items()}

            # 1-bp clipped versions (single-sided clipping only)
            seq_to_generated_header_1clip = {v[1:]: k for k, v in ref_generated_dict.items() if len(v) > 1}
            seq_to_provided_header_1clip = {v[1:]: k for k, v in ref_provided_dict.items() if len(v) > 1}

            generated_to_provided_headers = {}

            for seq in seqs_to_use:
                if seq in seq_to_generated_header and seq in seq_to_provided_header:
                    # Exact match
                    generated_to_provided_headers[seq_to_generated_header[seq]] = seq_to_provided_header[seq]
                elif seq in seq_to_provided_header_1clip:  
                    # Generated sequence has one extra base at the start
                    generated_to_provided_headers[seq_to_generated_header[seq]] = seq_to_provided_header_1clip[seq]
                elif seq in seq_to_generated_header_1clip:
                    # Provided sequence has one extra base at the start
                    generated_to_provided_headers[seq_to_generated_header_1clip[seq]] = seq_to_provided_header[seq]
                else:
                    print(f"Warning: No header mapping found for sequence: {seq}")

            # Now update index mappings
            panel_info_idx_to_provided_headers = {
                ref_generated_header_to_index[h]: generated_to_provided_headers[h] for h in generated_to_provided_headers
            }

            # Create markers table with corrected headers
            markers_table = createMarkersTable(panel_info_bed, panel_info_idx_to_provided_headers)

    else:
        # No comparison needed since reference.fasta is the auto-generated file.
        _, header_to_panel_info_index = createHeaderSeqDict("reference.fasta")
        panel_info_index_to_header = {v : k for k, v in header_to_panel_info_index.items()}
        markers_table = createMarkersTable(panel_info_bed, panel_info_index_to_header)
    
    # Check 2 - Check that fasta headers are similar
    ref_entries = readFasta(ref_file)
    fw_entries = readFasta(fw_file)
    rv_entries = readFasta(rv_file)

    cleaned_ref_headers = [cleanHeader(h) for h, _ in ref_entries]
    cleaned_fw_headers = [cleanHeader(h) for h, _ in fw_entries]
    cleaned_rv_headers = [cleanHeader(h) for h, _ in rv_entries]

    # Ensure sequence counts match
    if len(cleaned_ref_headers) != len(cleaned_fw_headers) or len(cleaned_ref_headers) != len(cleaned_rv_headers):
        raise ValueError("Mismatch in number of sequences across reference and primer files.")
    
    if cleaned_ref_headers == cleaned_fw_headers and cleaned_ref_headers == cleaned_rv_headers:
        print("Headers match across all files. No changes needed.")
    else:
        print("Headers mismatch detected. Overwriting primer files with corrected headers.")
        writeFasta(fw_file, cleaned_ref_headers, [seq for _, seq in fw_entries])
        writeFasta(rv_file, cleaned_ref_headers, [seq for _, seq in rv_entries])
        print(f"Headers corrected and overwritten in {fw_file} and {rv_file}.")

    # Create markers table
    print("Generating markers table...")
    markers_table.to_csv(markers_table_file, index=False)

def createHeaderSeqDict(fasta):
    index = 0
    index_dict = defaultdict()
    curr_fasta_dict = defaultdict(lambda : "")
    if os.path.isfile(fasta):
        headers_num_seen = defaultdict(lambda: 0) # Keep track of how many headers we have seen
        with open(fasta) as f:
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith(">"):
                    header = line[1:]
                    headers_num_seen[header] += 1
                    index_dict[header] = index
                    index += 1
                else:
                    if headers_num_seen[header] > 1:
                        raise Exception("Non-unique FASTA header detected. Every sequence must have a unique header.")
                    else:
                        curr_fasta_dict[header] += line 
    else:
       raise Exception(f"Non-unique FASTA header in {fasta} detected. Every sequence must have a unique header.")
    return curr_fasta_dict, index_dict

def createMarkersTable(bed, index_to_names):
    amplicon_bed = pd.read_csv(filepath_or_buffer=bed, sep=None, engine="python", header=None, names=["chromosome", "start", "end"])
    amplicon_bed['start'] += 1 #Need this for BED to Markers Table Conversion
    amplicon_bed['amplicon'] = [index_to_names[idx] for idx in amplicon_bed.index.tolist()]
    amplicon_bed['pos'] = (amplicon_bed['start'] + amplicon_bed['end'] ) // 2
    amplicon_bed['length'] = amplicon_bed['end'] - amplicon_bed['start'] + 1
    markers_table = amplicon_bed[['amplicon', 'chromosome', 'start', 'end', 'pos', 'length']]

    return markers_table

def cleanHeader(header):
    return re.split(r";()", header)[0]

def readFasta(file_path):
    entries = []
    with open(file_path, "r") as file:
        header, seq = None, []
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    entries.append((header, "".join(seq)))
                header, seq = line, []
            else:
                seq.append(line)
        if header:
            entries.append((header, "".join(seq)))  # Add last entry
    return entries

def writeFasta(file_path, headers, sequences):
    with open(file_path, "w") as file:
        for header, seq in zip(headers, sequences):
            file.write(f"{header}\n{seq}\n")

if __name__ == "__main__":
	main()
