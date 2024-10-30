"""
This module provides functionalities for reading the alignement and merging broken genes.
"""

####################################################################################################
# Imports
####################################################################################################

import os
import re
import numpy as np
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq

####################################################################################################
# Functions
####################################################################################################

def read_alignment(file_path, separator=";"):
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

def compare_sequences(sequence_one_seq, sequence_one_id, sequence_two, genome, output_folder, threshold_paralog=0.5, separator=";", force=False, multiple=False):
    """
    Compare two sequences and return the matching nucleotides.

    Args:
        sequence_one_seq (numpy.ndarray): The first sequence.
        sequence_one_id (str): The id of the first sequence.
        sequence_two (Bio.SeqRecord): The second sequence.
        genome (str): The name of the genome.
        output_folder (str): The path of the output folder.
        threshold_paralog (float): The threshold to consider the sequences as paralogs.
        separator (str): The separator to use to split the gene name.
        force (bool): Force the merge of the genes even if they might be paralogs.
        multiple (bool): If the sequences are multiple (more than 2 sequences).

    Returns:
        numpy.ndarray: An array containing the matching nucleotides.
    """

    # Convert the first sequence to a numpy array
    first_sequence = sequence_one_seq.copy()

    # Get the locations of the dashes and nucleotides in the first sequence
    dashes_locations = np.where(first_sequence == b"-")[0]
    nucleotides_locations = np.where(first_sequence != b"-")[0]

    tmp_sequence = np.frombuffer(str(sequence_two.seq).encode("utf-8"), dtype="S1")
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

    # Test for paralogs
    coverage_match_first_seq = position_matching.size / nucleotides_locations.size
    coverage_match_second_seq = position_matching.size / tmp_nucleotides_locations.size

    if coverage_match_first_seq >= threshold_paralog or coverage_match_second_seq >= threshold_paralog:

        with open(os.path.join(output_folder, "broken_genes.problems.paralogs.txt"), "a", encoding="utf-8") as file:
            sequence1 = b"".join(first_sequence).decode("utf-8")
            sequence2 = b"".join(tmp_sequence).decode("utf-8")

            # Find the matching nucleotides

            match_seq = b"".join(match_seq).decode("utf-8")

            message = f"Problem with {genome}, may be duplication:"
            subline = "="*len(message)

            name_seq1 = sequence_one_id
            name_seq2 = sequence_two.id

            length_name_seq1 = len(name_seq1)
            length_name_seq2 = len(name_seq2)

            if length_name_seq1 > length_name_seq2:
                name_seq2 += " "*(length_name_seq1 - length_name_seq2)
            elif length_name_seq2 > length_name_seq1:
                name_seq1 += " "*(length_name_seq2 - length_name_seq1)

            match_seq_name = " "*max(length_name_seq1, length_name_seq2)

            file.write("\n")
            file.write(f"{message}\n")
            file.write(f"{subline}\n\n")
            file.write(f"{name_seq1} {sequence1}\n")
            file.write(f"{match_seq_name} {match_seq}\n")
            file.write(f"{name_seq2} {sequence2}\n")

            return None

    if position_matching.size > 0:
        # Look if from the oposite side there is only gaps
        left_right = np.concatenate([first_sequence[:position_matching[0]], tmp_sequence[position_matching[-1]+1:]])
        right_left = np.concatenate([tmp_sequence[:position_matching[0]], first_sequence[position_matching[-1]+1:]])

        all_dashes_one_side = np.all(left_right == b'-') or np.all(right_left == b'-')
    else:
        all_dashes_one_side = False

    # Test if the location in front is empty (full of gaps)
    if np.all(tmp_sequence[nucleotides_locations] == b'-') and np.all(first_sequence[tmp_nucleotides_locations] == b'-'):
        first_sequence[dashes_locations] = tmp_sequence[dashes_locations]

        with open(os.path.join(output_folder, "broken_genes.merged.csv"), "a", encoding="utf-8") as file:
            if multiple:
                file.write(f"{genome}{separator}merged_multiple,{sequence_two.id} && {sequence_one_id}\n")
            else:
                file.write(f"{genome}{separator}merged,{sequence_one_id} && {sequence_two.id}\n")
    
    # Test if the matching position equal the matching length and if both the is gap after or before the matching position depending on the side
    elif (position_matching[-1] - position_matching[0] == len(position_matching) - 1) and (all_dashes_one_side):
        first_sequence[dashes_locations] = tmp_sequence[dashes_locations]

        with open(os.path.join(output_folder, "broken_genes.merged.csv"), "a", encoding="utf-8") as file:
            if multiple:
                file.write(f"{genome}{separator}merged_multiple,{sequence_one_id} && {sequence_two.id}\n")
            else:
                file.write(f"{genome}{separator}merged,{sequence_one_id} && {sequence_two.id}\n")
    
    # Everything else
    else:
        with open(os.path.join(output_folder, "broken_genes.problems.txt"), "a", encoding="utf-8") as file:
            sequence1 = b"".join(first_sequence).decode("utf-8")
            sequence2 = b"".join(tmp_sequence).decode("utf-8")

            # Find the matching nucleotides

            match_seq = b"".join(match_seq).decode("utf-8")

            message = f"Problem with {genome}, may be duplication:"
            subline = "="*len(message)

            name_seq1 = sequence_one_id
            name_seq2 = sequence_two.id

            length_name_seq1 = len(name_seq1)
            length_name_seq2 = len(name_seq2)

            if length_name_seq1 > length_name_seq2:
                name_seq2 += " "*(length_name_seq1 - length_name_seq2)
            elif length_name_seq2 > length_name_seq1:
                name_seq1 += " "*(length_name_seq2 - length_name_seq1)

            match_seq_name = " "*max(length_name_seq1, length_name_seq2)

            file.write("\n")
            file.write(f"{message}\n")
            file.write(f"{subline}\n\n")
            file.write(f"{name_seq1} {sequence1}\n")
            file.write(f"{match_seq_name} {match_seq}\n")
            file.write(f"{name_seq2} {sequence2}\n")

            if force:
                # Test if the second sequence is before the first sequence else the opposite
                if tmp_nucleotides_locations[0] < nucleotides_locations[0]:
                    first_sequence = np.concatenate([tmp_sequence[tmp_nucleotides_locations], first_sequence[nucleotides_locations[-1]+1:]])

                    if not multiple:
                        with open(os.path.join(output_folder, "broken_genes.merged.csv"), "a", encoding="utf-8") as file:
                            file.write(f"{genome}{separator}merged,{sequence_two.id} && {sequence_one_id}\n")
                else:
                    first_sequence = np.concatenate([first_sequence[nucleotides_locations], tmp_sequence[tmp_nucleotides_locations[-1]+1:]])

                    if not multiple:
                        with open(os.path.join(output_folder, "broken_genes.merged.csv"), "a", encoding="utf-8") as file:
                            file.write(f"{genome}{separator}merged,{sequence_one_id} && {sequence_two.id}\n")

    return first_sequence

####################################################################################################

def merge_broken_genes(alignment_dict, output_folder, separator=";", force=False):
    """
    Merge broken genes in an alignment.

    Args:
        alignment (dict): A dict of SeqRecord objects representing the sequences in the alignment.
        output_folder (str): The path of the output folder.
        separator (str): The separator to use to split the gene name.
        force (bool): Force the merge of the genes even if they might be paralogs.
        
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

        if len(value) == 2:
            # Get the the first sequence as pivot
            sequence_one = np.frombuffer(str(value[0].seq).encode("utf-8"), dtype="S1").copy()
            sequence_two = value[1]

            concatenate_sequence = compare_sequences(
                sequence_one_seq=sequence_one,
                sequence_one_id=value[0].id,
                sequence_two=sequence_two,
                genome=key,
                output_folder=output_folder,
                separator=separator,
                force=force,
                multiple=False,
            )

        # More than 2 sequences so complicated case
        else:
            # Sort the sequences based on the start of the sequences
            sorted_sequences = sorted(value, key=lambda x: re.search(r'[^-]', str(x.seq)).start())

            # Get the the first sequence as pivot
            concatenate_sequence = np.frombuffer(str(sorted_sequences[0].seq).encode("utf-8"), dtype="S1").copy()
            concatenate_sequence_id = sorted_sequences[0].id
            first_sequence_len = len(sorted_sequences[0].seq)

            for sequence_two in sorted_sequences[1:]:
                concatenate_sequence = compare_sequences(
                    sequence_one_seq=concatenate_sequence,
                    sequence_one_id=concatenate_sequence_id,
                    sequence_two=sequence_two,
                    genome=key,
                    output_folder=output_folder,
                    separator=separator,
                    force=force,
                    multiple=True,
                )
                
                # Because the sequence is now a merged sequence
                concatenate_sequence_id = f"{key}{separator}merged_multiple"

                # Break the loop if paralog
                if concatenate_sequence is None:
                    break
                # Break you had to concatenate two sequence one after the other
                elif concatenate_sequence.size != first_sequence_len:
                    concatenate_sequence = [np.frombuffer(str(sequence.seq).encode("utf-8"), dtype="S1").copy() for sequence in sorted_sequences]

                    with open(os.path.join(output_folder, "broken_genes.merged.csv"), "a", encoding="utf-8") as file:
                        file.write(f"{key}{separator}merged_multiple,{' && '.join([sequence.id for sequence in sorted_sequences])}\n")
                    break

        # Check if the sequence is not empty
        if concatenate_sequence is not None:
            # Update the sequence in the alignment
            tmp_record = value[0]
            tmp_record.seq = Seq(b"".join(concatenate_sequence).decode("utf-8"))
            tmp_record.id = f"{key}{separator}merged"
            tmp_record.description = ""
            tmp_record.name = ""

            new_alignment[key] = [tmp_record]

    alignment_length_1.update(new_alignment)

    return alignment_length_1

####################################################################################################

def write_alignment(alignment, output_file, force=False):
    """
    Write the alignment to a file.

    Args:
        alignment (dict): A dict of SeqRecord objects representing the sequences in the alignment.
        output_file (str): The path of the output file.
        force (bool): Force the merge of the genes even if they might be paral
    """
    with open(output_file, "w", encoding="utf-8") as handle:
        for key, value in alignment.items():
            for record in value:
                if force:
                    record.seq = record.seq.replace("-", "")

                SeqIO.write(record, handle, "fasta")

####################################################################################################