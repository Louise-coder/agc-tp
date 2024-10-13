#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import gzip
import textwrap
from pathlib import Path
from typing import Iterator, List

# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Louise LAM"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Louise LAM"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Louise LAM"
__email__ = "louise.lam@etu.u-paris.fr"
__status__ = "Developpement"


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments():  # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(
        description=__doc__, usage=f"{sys.argv[0]} -h"
    )
    parser.add_argument(
        "-i",
        "-amplicon_file",
        dest="amplicon_file",
        type=isfile,
        required=True,
        help="Amplicon is a compressed fasta file (.fasta.gz)",
    )
    parser.add_argument(
        "-s",
        "-minseqlen",
        dest="minseqlen",
        type=int,
        default=400,
        help="Minimum sequence length for dereplication (default 400)",
    )
    parser.add_argument(
        "-m",
        "-mincount",
        dest="mincount",
        type=int,
        default=10,
        help="Minimum count for dereplication  (default 10)",
    )
    parser.add_argument(
        "-o",
        "-output_file",
        dest="output_file",
        type=Path,
        default=Path("OTU.fasta"),
        help="Output file",
    )
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """
    seq = ""
    with gzip.open(amplicon_file, "rt") as monfich:
        for line in monfich:
            if line.startswith(">"):
                if len(seq) >= minseqlen:
                    yield seq
                seq = ""
            else:
                seq += line.strip()
        if len(seq) >= minseqlen:  # last sequence
            yield seq


def dereplication_fulllength(
    amplicon_file: Path, minseqlen: int, mincount: int
) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of sequence
    with a count >= mincount and a length >= minseqlen.
    """
    counts = {}
    for seq in read_fasta(amplicon_file, minseqlen):
        counts[seq] = counts.get(seq, 0) + 1
    sorted_counts = sorted(
        counts.items(), key=lambda item: item[1], reverse=True
    )
    for seq, count in sorted_counts:
        if count >= mincount:
            yield [seq, count]


def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences in the
    format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    seq_a = alignment_list[0]
    seq_b = alignment_list[1]
    identity = 0
    for i in range(len(seq_a)):
        if seq_a[i] == seq_b[i]:
            identity += 1
    return identity * 100 / len(seq_a)


def abundance_greedy_clustering(
    amplicon_file: Path,
    minseqlen: int,
    mincount: int,
    chunk_size: int,
    kmer_size: int,
) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    most_abundant = dereplication_fulllength(
        amplicon_file, minseqlen, mincount
    )
    all_otus = []
    for seq, count in most_abundant:
        is_otu = all(
            get_identity(
                nw.global_align(
                    seq,
                    otu_seq,
                    gap_open=-1,
                    gap_extend=-1,
                    matrix=str(Path(__file__).parent / "MATCH"),
                )
            )
            <= 97
            for otu_seq, _ in all_otus
        )
        if is_otu:
            all_otus.append([seq, count])
    return all_otus


def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, "w", encoding="utf-8") as file:
        for i, (otu, count) in enumerate(OTU_list):
            file.write(f">OTU_{i+1} occurrence:{count}\n")
            file.write(f"{textwrap.fill(otu, 80)}\n")


# ==============================================================
# Main program
# ==============================================================
def main():  # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici
    all_otus = abundance_greedy_clustering(
        args.amplicon_file, args.minseqlen, args.mincount, 0, 0
    )
    write_OTU(all_otus, args.output_file)


if __name__ == "__main__":
    main()
