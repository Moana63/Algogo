from collections import Counter
from argparse import ArgumentParser
from os import getpid
from psutil import Process
from random import choice
from itertools import product


def frequency(read: str, seed_size) -> Counter:
    """Returns the frequency of kmers per read
    Parameters
    ----------
    read : str
        a DNA read
    seed_size : int
        size of the kmer
    Returns
    -------
    dict
        a dictionnary containing the kmers encountered in the sequence and their number of occurences
    """
    return Counter([read[k:k+seed_size] for k in range(len(read)-seed_size)])


def binary(dico: dict, len_read: int, list_xmers: list, threshold: float) -> str:
    """Take the proportions of each possible kmers in the sequence and returns a sequence of 0 and 1. 
    For each kmer proportion, if the proportion is superior to the threshold, a 1 is added to the sequence. Else, a 0 is added.

    Parameters
    ----------
    dico : dict
        a dictionnary containing the proportions of each different kmers encountered in a sequence.
    len_read : int
        The length of the sequence
    list_xmers : list
        the list of all possible kmers of size k
    threshold : float
        the threshold above which the kmer proportion will be put to 1 and under which it will be put to 0.
    Returns
    -------
    str
        the sequence of 0 and 1 as a simplification of the proportions of kmers
    """
    binary_seq = str()
    for kmer in list_xmers:
        if kmer in dico.keys():
            if dico[kmer]/(len_read-2) <= threshold:
                binary_seq += "0"
            else:
                binary_seq += "1"
        else:
            binary_seq += "0"
    return binary_seq


def indexation(list_seq: list, seed_size: int, len_read: int = 100) -> dict:
    """Index all the read sequences from the fasta file according to their identifier. 
    The identifier is a sequence of 0 and 1 linked to the proportions of différent kmers in the read sequence.

    Parameters
    ----------
    list_seq : list
        a list containing all the read sequence from the fasta file
    seed_size : int
        the size k of the kmer
    len_read : int, optional
        the length of the read sequences in the file, by default 100

    Returns
    -------
    dict
        an index with the identifier sequence of 0 and 1 as key and a list of the indexes in list_seq of the corresponding reads as values. 
    """
    index = dict()
    list_xmers = ["".join(tuples)
                  for tuples in list(product('ATCG', repeat=seed_size))]
    # The threshold is calculated considering that the occurence of each kmer is equiprobable
    # If 4 differents kmers exist, then each kmer have a propability of 25% percent to occur. The threshold will be set to 0.25.
    threshold = 1/len(list_xmers)
    for i, read in enumerate(list_seq):
        dico = frequency(read, seed_size)
        binary_seq = binary(dico, len_read, list_xmers, threshold)
        if binary_seq in index:
            index[binary_seq] += [i]
        else:
            index[binary_seq] = [i]
    return index


def write_file(filein, fileout, seed_size: int = 4):
    list_seq = clean_fasta(filein)
    index = indexation(list_seq, seed_size)
    list_keys = sorted(index.keys())
    with open(fileout, "w") as out_file:
        for key in list_keys:
            for seq in index[key]:
                out_file.write(list_seq[int(seq)]+'\n')
