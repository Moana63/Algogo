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


def binary(dico, len_read, list_xmers, seuil):
    binary_seq = str()
    for kmer in list_xmers:
        if kmer in dico.keys():
            if dico[kmer]/(len_read-2) <= seuil:
                binary_seq += "0"
            else:
                binary_seq += "1"
        else:
            binary_seq += "0"
    return binary_seq


def indexation(list_seq, seed_size, len_read=100):
    index = dict()
    list_xmers = ["".join(tuples)
                  for tuples in list(product('ATCG', repeat=seed_size))]
    seuil = 1/len(list_xmers)
    for i, read in enumerate(list_seq):
        dico = frequency(read, seed_size)
        binary_seq = binary(dico, len_read, list_xmers, seuil)
        if binary_seq in index:
            index[binary_seq] += [i]
        else:
            index[binary_seq] = [i]
    return index


def write_file(filein, fileout, seed_size=2):
    with open(filein) as f:
        list_seq = [seq.strip() for seq in f]
    index = indexation(list_seq, seed_size)
    list_keys = sorted(index.keys())
    with open(fileout, "w") as out_file:
        for key in list_keys:
            for seq in index[key]:
                out_file.write(list_seq[int(seq)]+'\n')


# write_file("clean_file_test", "clean_file_out_proba_2mers")
