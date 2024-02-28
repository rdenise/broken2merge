"""
This module provides functionalities for reading the alignement and merging broken genes.
"""

####################################################################################################
# Imports
####################################################################################################

import os
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from tqdm import tqdm

####################################################################################################
# Functions
####################################################################################################

def read_alignment(file_path, separator="_"):
    """
    Read an alignment in fasta format using Biopython.

    Args:
        file_path (str): The path to the alignment file.
        separator (str): The separator to use to split the gene name.

    Returns:
        list: A list of SeqRecord objects representing the sequences in the alignment.
    """
    alignment = {}

    with open(file_path, "r", encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            record_id = record.id.split(separator)[0]

            if record_id not in alignment:
                alignment[record_id] = [record]
            else:
                alignment[record_id].append(record)

    return alignment

####################################################################################################

def split_dictionary(dictionary):
    """
    Split a dictionary into two dictionaries based on the length of the values.

    Args:
        dictionary (dict): The dictionary to be split.

    Returns:
        tuple: A tuple containing two dictionaries. The first dictionary contains the key-value pairs
               where the length of the value is 1, and the second dictionary contains the key-value pairs
               where the length of the value is greater than 1.
    """
    dictionary_length_1 = {}
    dictionary_length_greater_than_1 = {}

    for key, value in dictionary.items():
        if len(value) == 1:
            dictionary_length_1[key] = value
        else:
            dictionary_length_greater_than_1[key] = value

    return dictionary_length_1, dictionary_length_greater_than_1

####################################################################################################

def merge_broken_genes(alignment_dict, output_folder, separator="_"):
    """
    Merge broken genes in an alignment.

    Args:
        alignment (dict): A dict of SeqRecord objects representing the sequences in the alignment.
        output_folder (str): The path of the output folder.
        separator (str): The separator to use to split the gene name.
        
    Returns:
        list: A dict of SeqRecord objects representing the sequences in the alignment after merging
              the broken genes.
    """

    # Split the alignment into two dictionaries based on the length of the values
    alignment_length_1, alignment_length_greater_than_1 = split_dictionary(alignment_dict)

    with open(os.path.join(output_folder, "broken_genes.problems.txt"), "w", encoding="utf-8") as file:
        file.write("**** Please check for problems ****\n")

    with open(os.path.join(output_folder, "broken_genes.merged.csv"), "w", encoding="utf-8") as file:
        file.write("Merged_seq_name,seqid_used\n")

    new_alignment = {}

    # Merge the broken genes
    for key, value in tqdm(alignment_length_greater_than_1.items(), desc="Parsing duplicate names"):
        # Get the the first sequence as pivot
        first_sequence = value[0]

        # Convert the first sequence to a numpy array
        first_sequence = np.frombuffer(str(first_sequence.seq).encode("utf-8"), dtype="S1").copy()

        # Get the locations of the dashes and nucleotides in the first sequence
        dashes_locations = np.where(first_sequence == b"-")[0]
        nucleotides_locations = np.where(first_sequence != b"-")[0]

        # Iterate over the other sequences and merge the broken genes
        for sequence in value[1:]:
            tmp_sequence = np.frombuffer(str(sequence.seq).encode("utf-8"), dtype="S1")
            tmp_dash_locations = np.where(tmp_sequence == b"-")[0]
            tmp_nucleotides_locations = np.where(tmp_sequence != b"-")[0]

            # print('tmp_sequence')
            # print(b"".join(tmp_sequence[nucleotides_locations]).decode("utf-8"))
            # print(np.all(tmp_sequence[nucleotides_locations] == b'-'))
            # print('first_sequence')
            # print(b"".join(first_sequence[tmp_nucleotides_locations]).decode("utf-8"))
            # print(np.all(first_sequence[tmp_nucleotides_locations] == b'-'))
            
            # Find the matching nucleotides and removing matching gaps
            matching = np.where(first_sequence == tmp_sequence)[0]
            match_seq = np.array([b' ']*first_sequence.shape[0], dtype='|S1')
            match_seq[matching] = b'|'
            match_seq[dashes_locations] = b' '

            # Get the positions of the matching nucleotides (np.where return index)
            position_matching = np.where(match_seq == b'|')[0]

            if position_matching.size > 0:
                # Look if from the oposite side there is only gaps
                left_right = np.concatenate([first_sequence[:position_matching[0]], tmp_sequence[position_matching[-1]+1:]])
                right_left = np.concatenate([tmp_sequence[:position_matching[0]], first_sequence[position_matching[-1]+1:]])

                all_dashes_one_side = np.all(left_right == b'-') or np.all(right_left == b'-')
            else:
                all_dashes_one_side = False

            if np.all(tmp_sequence[nucleotides_locations] == b'-') and np.all(first_sequence[tmp_nucleotides_locations] == b'-'):
                first_sequence[dashes_locations] = tmp_sequence[dashes_locations]

                with open(os.path.join(output_folder, "broken_genes.merged.csv"), "a", encoding="utf-8") as file:
                    file.write(f"{key}{separator}merged,{value[0].id} && {sequence.id}\n")
            
            elif (position_matching[-1] - position_matching[0] == len(position_matching) - 1) and (all_dashes_one_side):
                first_sequence[dashes_locations] = tmp_sequence[dashes_locations]

                with open(os.path.join(output_folder, "broken_genes.merged.csv"), "a", encoding="utf-8") as file:
                    file.write(f"{key}{separator}merged,{value[0].id} && {sequence.id}\n")
            else:
                with open(os.path.join(output_folder, "broken_genes.problems.txt"), "a", encoding="utf-8") as file:
                    sequence1 = b"".join(first_sequence).decode("utf-8")
                    sequence2 = b"".join(tmp_sequence).decode("utf-8")

                    # Find the matching nucleotides

                    match_seq = b"".join(match_seq).decode("utf-8")

                    message = f"Problem with {key}, may be duplication:"
                    subline = "="*len(message)

                    file.write("\n")
                    file.write(f"{message}\n")
                    file.write(f"{subline}\n\n")
                    file.write(f"{sequence1}\n")
                    file.write(f"{match_seq}\n")
                    file.write(f"{sequence2}\n")

        # Update the sequence in the alignment
        tmp_record = value[0]
        tmp_record.seq = Seq(b"".join(first_sequence).decode("utf-8"))
        tmp_record.id = f"{key}{separator}merged"
        tmp_record.description = ""
        tmp_record.name = ""

        new_alignment[key] = [tmp_record]

    alignment_length_1.update(new_alignment)

    return alignment_length_1

####################################################################################################

def write_alignment(alignment, output_file):
    """
    Write the alignment to a file.

    Args:
        alignment (dict): A dict of SeqRecord objects representing the sequences in the alignment.
        output_file (str): The path of the output file.
    """
    with open(output_file, "w", encoding="utf-8") as handle:
        for key, value in alignment.items():
            for record in value:
                SeqIO.write(record, handle, "fasta")

####################################################################################################