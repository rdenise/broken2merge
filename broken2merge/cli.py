# -*- coding: utf-8 -*-

"""
This module is used to handle command line interface (CLI) for the taxmyphage package.
It includes classes and functions to parse command line arguments and display help messages.
"""

####################################################################################################
# Imports
####################################################################################################

import os
from argparse import ArgumentParser
from broken2merge import __version__

####################################################################################################
# Classes
####################################################################################################

def cli(args=None):
    """
    Command line interface for taxmyphage
    """

    description = """Takes a alignement and find genes that seems broken and merge them together"""


    parser = ArgumentParser(
        description=description, conflict_handler="resolve", prog="broken_merge"
    )

    ####################################################################################################
    # Create general subparser that will be given to all subparsers
    ####################################################################################################

    general_options = parser.add_argument_group(title="General input dataset options")

    general_options.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Show verbose output. (For debugging purposes)",
    )
 
    general_options.add_argument(
        "-V",
        "--version",
        action="version",
        help="Show the version number and exit.",
        version=f"broken_merge v{__version__}",
    )

    general_options.add_argument(
        "-i",
        "--input",
        dest="in_fasta",
        metavar="FASTA_FILE",
        type=str,
        help="Path to an input fasta file (Required)",
        required=True,
    )

    general_options.add_argument(
        "-o",
        "--output",
        type=str,
        default=os.path.join(os.getcwd(), "merge_broken_res"),
        dest="output",
        help="Path of the output folder (Default: merge_broken_res)",
    )

    general_options.add_argument(
        "-s",
        "--separator",
        type=str,
        default=";",
        help="Separator to use to split the gene name (Default: ';')",
    )

    general_options.add_argument(
        "--force_merge",
        default=False,
        action="store_true",
        help="Force the merge of the genes even if they might be paralogs, will give a unaligned file as main input",
    )

    # general_options.add_argument(
    #     "--discard",
    #     action="store_true",
    #     type=bool,
    #     default=False,
    #     help="Discard the sequence from the alignment in case of duplication",
    # )

    args, nargs = parser.parse_known_args(args)

    return args, nargs
    # No return value means no error.
    # Return a value of 1 or higher to signify an error.
    # See https://docs.python.org/3/library/sys.html#sys.exit


####################################################################################################
