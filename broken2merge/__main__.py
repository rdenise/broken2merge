#!/usr/bin/env python3

"""
This module is used for assembling into one genes that seems broken
"""

###################################################################################################
# Imports
###################################################################################################


import os

from broken2merge import cli
from broken2merge.utils import create_folder
from broken2merge.alignment import merge_broken_genes, read_alignment, write_alignment

###################################################################################################
# Functions
###################################################################################################


def main():
    """
    This is the main function of the script
    """
    # Set up the arguments
    args, nargs = cli.cli()

    verbose = args.verbose

    create_folder(args.output)

    file_extension = os.path.splitext(args.in_fasta)[1]

    dict_alignment = read_alignment(file_path=args.in_fasta,
                                    separator=args.separator)

    merged_alignment = merge_broken_genes(alignment_dict=dict_alignment,
                                          output_folder=args.output,
                                          separator=args.separator,
                                          force=args.force_merge,
                                          )

    basename_alignment = os.path.basename(args.in_fasta)

    final_alignment = os.path.join(args.output, basename_alignment.replace(file_extension, ".merged" + file_extension))

    write_alignment(alignment=merged_alignment, output_file=final_alignment, force=args.force_merge)

###################################################################################################

if __name__ == "__main__":
    main()

###################################################################################################
