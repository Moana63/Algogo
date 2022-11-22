from collections import Counter
from argparse import ArgumentParser
from os import getpid
from psutil import Process
from random import choice


def frequency(read: str, seed_size=2) -> Counter:
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


def binary(dico, len_read, list_2mers, seuil):
    binary_seq = str()
    for kmer in list_2mers:
        if kmer in dico.keys():
            if dico[kmer]/(len_read-2) <= seuil:
                binary_seq += "0"
            else:
                binary_seq += "1"
        else:
            binary_seq += "0"
    return binary_seq


def indexation(list_seq, len_read=100):
    index = dict()
    list_2mers = ["".join([i, j]) for i in ['A', 'T', 'C', 'G']
                  for j in ['A', 'T', 'C', 'G']]
    seuil = 1/len(list_2mers)
    for i, read in enumerate(list_seq):
        dico = frequency(read, 2)
        binary_seq = binary(dico, len_read, list_2mers, seuil)
        if binary_seq in index:
            index[binary_seq] += [i]
        else:
            index[binary_seq] = [i]
    return index


def write_file(filein, fileout):
    with open(filein) as f:
        list_seq = [seq.strip() for seq in f]
    index = indexation(list_seq)
    with open(fileout, "w") as out_file:
        for group in index.values():
            for seq in group:
                out_file.write(list_seq[int(seq)]+'\n')


write_file("clean_file_test", "clean_file_out_proba_2mers")
