from collections import Counter
from argparse import ArgumentParser
from os import getpid
from psutil import Process
from random import choice
from itertools import product


def frequency_minimizer(read: str, seed_size: int, len_window: int = 10) -> dict:
    """Returns the frequency of minimisers per read

    Parameters
    ----------
    read : str
        a DNA read
    seed_size : int
        size of the minimiser
    len_window : int, optional
        length of the window sliding on the sequence, by default 10

    Returns
    -------
    dict
        a dictionnary containing the minimisers encountered in the sequence as key and their number of occurences as values
    """
    list_minimisers = list()
    i = 0
    while i < (len(read)-len_window+1):
        minimiser = min([(read[i+j: i+j+seed_size], j)
                        for j in range(len_window-seed_size+1)])
        list_minimisers.append(minimiser[0])
        # we jump to the next window not containing the last minimiser encountered, to reduce computation time
        i += minimiser[1] + 1
    return Counter(list_minimisers)


def binary_minimisers(dico: dict, list_xmers: list) -> str:
    """Take the proportions of each possible minimisers in the sequence and returns a sequence of 0 and 1. 
    For each minimisers proportion, if the minimiser is present, a 1 is added to the sequence. Else, a 0 is added.

    Parameters
    ----------
    dico : dict
        a dictionnary containing the proportions of each different minimisers encountered in a sequence.
    list_xmers : list
        the list of all possible minimisers

    Returns
    -------
    str
        the sequence of 0 and 1 as a simplification of the presence and absence of minimisers
    """

    binary_seq = str()
    for kmer in list_xmers:
        if kmer in dico.keys():
            binary_seq += "1"
        else:
            binary_seq += "0"
    return binary_seq


def indexation_minimisers(list_seq: list, seed_size: int) -> dict:
    """Index all the read sequences from the fasta file according to their identifier. 
    The identifier is a sequence of 0 and 1 linked to the presence of different minimisers in the read sequence.

    Parameters
    ----------
    list_seq : list
        a list containing all the read sequence from the fasta file
    seed_size : int
       the size of the minimiser

    Returns
    -------
    dict
        an index with the identifier sequence of 0 and 1 as key and a list of the indexes in list_seq of the corresponding reads as values. 
    """
    index = dict()
    list_xmers = sorted(["".join(tuples)
                         for tuples in list(product('ATCG', repeat=seed_size))])
    for i, read in enumerate(list_seq):
        dico = frequency_minimizer(read, seed_size)
        binary_seq = binary_minimisers(dico, list_xmers)
        if binary_seq in index:
            index[binary_seq] += [i]
        else:
            index[binary_seq] = [i]
    return index


def write_file(filein, fileout, seed_size: int = 4, len_window: int = 33):
    list_seq = clean_fasta(filein)
    index = indexation_minimisers(list_seq, seed_size, len_window)
    list_keys = sorted(index.keys())
    with open(fileout, "w") as out_file:
        for key in list_keys:
            for seq in index[key]:
                out_file.write(list_seq[int(seq)]+'\n')
