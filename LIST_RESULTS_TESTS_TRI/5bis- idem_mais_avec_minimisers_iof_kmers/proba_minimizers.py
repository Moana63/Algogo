from collections import Counter
from argparse import ArgumentParser
from os import getpid
from psutil import Process
from random import choice
from itertools import product


def frequency(read: str, seed_size) -> Counter:
    """Returns kmer frequency per read

    Args:
        read (str): a DNA lecture
        seed_size (int, optional): Size of window. Defaults to 2.

    Returns:
        Counter: counts of kmers inside the lecture
    """
    return Counter([read[k:k+seed_size] for k in range(len(read)-seed_size)])


file = ["GCAACTCTGCTACGTCTGGGGCTGGCTTACGGCCCCGGGGGGATGTCATTACGTGAAGTCACTGCATGGGCTCAGCTCCATGACGTTGCAACATTATCTG",
        "ATGTGACAAGTTACGACGCCAGAGCCAGCACCTCAGCCACCCACTGCTGGAAGCCTTTCACATCCAGATCCAATGCCACCTGTACATTGGCTGGCTTGCC", "CGAATTCGACTACTGTTGCGTACACGCCTCGCTGGCGCTGCGCGAAGACGGTTACGAAACCATTATGGTTAACTGTAACCCGGAAACCGTCTCCACCGAC", "CGAATTCGACTACTGTTGCGTACACGCCTCGCTGGCGCTGCGCGAAGACGGTTACGAAACCATTATGGTTAACTGTCCCGGAAACCGTCTCCACCGACAA"]

# même principa que la fréquence de kmer, mais "filtrée" en ne comptant que la fréqence des minimsers
# 0 si minimiser non présent, 1 si présent


def frequency_minimizer(read, seed_size, len_window=10):
    list_minimisers = list()
    i = 0
    while i < (len(read)-len_window+1):
        minimiser = min([(read[i+j: i+j+seed_size], j)
                        for j in range(len_window-seed_size+1)])
        list_minimisers.append(minimiser[0])
        # on saute les positions de i pour dépasser le minimiser
        i += minimiser[1] + 1
    return Counter(list_minimisers)


def binary_minimisers(dico, list_xmers):
    binary_seq = str()
    for kmer in list_xmers:
        if kmer in dico.keys():
            binary_seq += "1"
        else:
            binary_seq += "0"
    return binary_seq


def indexation_minimisers(list_seq, seed_size):
    index = dict()
    list_xmers = ["".join(tuples)
                  for tuples in list(product('ATCG', repeat=seed_size))]
    for i, read in enumerate(list_seq):
        dico = frequency_minimizer(read, seed_size)
        binary_seq = binary_minimisers(dico, list_xmers)
        if binary_seq in index:
            index[binary_seq] += [i]
        else:
            index[binary_seq] = [i]
    return index


def write_file(filein, fileout, seed_size=2):
    with open(filein) as f:
        list_seq = [seq.strip() for seq in f]
    index = indexation_minimisers(list_seq, seed_size)
    list_keys = sorted(index.keys())
    with open(fileout, "w") as out_file:
        for key in list_keys:
            for seq in index[key]:
                out_file.write(list_seq[int(seq)]+'\n')


#write_file("clean_file", "clean_file_out_proba_2mers")
