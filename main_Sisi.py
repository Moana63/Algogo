from argparse import ArgumentParser
from os import getpid
from psutil import Process
from random import choice
from collections import Counter
from itertools import product
from operator import itemgetter

def write_output(output:str,reads_ordered:list) -> None:
    with open(output, 'w') as writer:
        writer.write('\n'.join(reads_ordered))

def frequency(read: str, seed_size=2) -> Counter:
    """Returns kmer frequency per read

    Args:
        read (str): a DNA lecture
        seed_size (int, optional): Size of window. Defaults to 2.

    Returns:
        Counter: counts of kmers inside the lecture
    """
    return Counter([read[k:k+seed_size] for k in range(len(read)-seed_size)])

######## KMER FREQUENCY ########

def compress_by_kmer(input: str, output: str, ksize: int = 4, kmer_number:int = 3) -> float:
    """Sort a read file by the most present kmers each read contains

    Args:
        input (str): input file
        output (str): output file
        ksize (int, optional): size of kmer. Defaults to 4.
        kmer_number (int, optional): number of top common kmers. Defaults to 3.

    Returns:
        float: _description_
    """
    init_memory = get_memory()
    reads: list = clean_fasta(input)
    vectors: list = [''.join([key for key,value in frequency(read, ksize).most_common(kmer_number)]) for read in reads]
    ordered_reads: list = [reads[i] for i in [i for i,_ in sorted(enumerate(vectors), key=lambda x:x[1])]]
    max_memory = get_memory()
    write_output(output,ordered_reads)
    return abs(max_memory - init_memory)

######## DISTANCE MATRIX => O(n**1000000000000) environ ########

def hash_kmer(kmer: str) -> int:
    return sum(({'A': 1, 'C': 2, 'G': 3, 'T': 4, }.get(l, 0)*(5**i) for i, l in enumerate(kmer)))


def calculate_algebric_distance(a: Counter, b: Counter) -> Counter:
    return sum([abs(a[elt]-b[elt]) for elt in set().union(*a, *b)])


def algebric_clustering(input: str, output: str, ksize: int = 4):
    init_memory = get_memory()
    reads: list = clean_fasta(input)
    ref: Counter = frequency(reads[0], ksize)
    reads_ordered: list = sorted(
        reads, key=lambda read: calculate_algebric_distance(frequency(read, ksize), ref))
    max_memory = get_memory()
    write_output(output,reads_ordered)
    return abs(max_memory - init_memory)


def insert(list_to_return: list, element_to_insert: int, position_to_insert: int) -> list:
    return [*list_to_return[:position_to_insert], element_to_insert, *list_to_return[position_to_insert:]]


def find_insert_point(distances: list):
    distances = [0, *distances, 0]
    pos_value: int = float('inf')
    pos: int = 0
    for i in range(len(distances)-1):
        if pos_value > distances[i] + distances[i+1]:
            pos = i+1
    return pos


def matrix_clustering(input: str, output: str, ksize: int = 4):
    init_memory = get_memory()
    reads: list = clean_fasta(input)
    vectors: list = [frequency(read, ksize) for read in reads]
    distance_matrix: list = [[calculate_algebric_distance(
        vectors[i], vectors[j]) if i != j else float('inf') for i in range(len(vectors))] for j in range(len(vectors))]
    sorted_reads: list = [0]  # we init with the first read from the list

    for i, distance in enumerate(distance_matrix):
        sorted_reads = insert(sorted_reads, i, find_insert_point(
            [distance[sorted_index] for sorted_index in sorted_reads]))

    ordered_reads: list = [reads[i] for i in sorted_reads]

    max_memory = get_memory()
    write_output(output,ordered_reads)
    return abs(max_memory - init_memory)


def read_to_vector(ref_vector, read, ksize) -> list:
    return []


def matrix_clustering_new(input: str, output: str, ksize: int = 4):
    ref_vector: list = [code for code in map(
        ''.join, product('ATCG', repeat=ksize))]
    init_memory = get_memory()
    reads: list = clean_fasta(input)
    vectors: list = [frequency(read, ksize) for read in reads]
    distance_matrix: list = [[calculate_algebric_distance(
        vectors[i], vectors[j]) if i > j else float('inf') for i in range(len(vectors))] for j in range(len(vectors))]
    sorted_reads: list = [0]  # we init with the first read from the list

    for i, distance in enumerate(distance_matrix):
        sorted_reads = insert(sorted_reads, i, find_insert_point(
            [distance[sorted_index] for sorted_index in sorted_reads]))

    ordered_reads: list = [reads[i] for i in sorted_reads]

    max_memory = get_memory()
    write_output(output,ordered_reads)
    return abs(max_memory - init_memory)


def get_memory() -> float:
    return Process(getpid()).memory_info().rss / 1024 ** 2


def clean_fasta(input: str) -> list[str]:
    with open(input, 'r') as reader:
        return [l.strip() for l in reader if l[0] not in ['\n', '>']]


def compress_naive(input: str, output: str):
    init_memory = get_memory()
    lines = sorted(clean_fasta(input))
    max_memory = get_memory()
    write_output(output,lines)
    return abs(max_memory - init_memory)


if __name__ == "__main__":
    """
    parser = ArgumentParser()
    parser.add_argument("-i", "--input",
                        help="Path to a fasta-like file to re-order", type=str)
    parser.add_argument("-o", "--output",
                        help="Path to a txt file, output of the program", type=str)
    args = parser.parse_args()
    """
    """
    kmers = [''.join(choice(['A', 'T', 'C', 'G'])
                     for _ in range(100)) for _ in range(1000000)]
    hashes = {k: hash_kmer(k) for k in kmers}
    print(len(hashes.keys()) == len(set(l for l in hashes.values())))
    """
    #compress_by_kmer("Ecoli_100Kb/ecoli_100Kb_reads_005x.fasta","",4)
