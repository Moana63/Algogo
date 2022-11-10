from argparse import ArgumentParser
from os import getpid
from psutil import Process
from random import choice
from collections import Counter




def hash_kmer(kmer: str) -> int:
    return sum(({'A': 1, 'C': 2, 'G': 3, 'T': 4, }.get(l, 0)*(5**i) for i, l in enumerate(kmer)))

def frequency(read:str,seed_size) -> Counter:
    return Counter([read[k:k+seed_size] for k in range(len(read)-seed_size)])

def calculate_algebric_distance(a:Counter,b:Counter) -> Counter:
    return sum([abs(a[elt]-b[elt]) for elt in set().union(*a,*b)])

def algebric_clustering(input: str, output: str, ksize:int=4):
    init_memory = get_memory()
    reads:list = clean_fasta(input)
    ref:Counter = frequency(reads[0],ksize)
    reads_ordered:list = sorted(reads, key=lambda read: calculate_algebric_distance(frequency(read,ksize), ref))
    max_memory = get_memory()
    with open(output, 'w') as writer:
        writer.write(''.join(reads_ordered))
    return abs(max_memory - init_memory)
    #vectors:list = [frequency(read,ksize) for read in reads]
    #while(sum(is_collected) != len(is_collected)):
    #distance_matrix:list = [[calculate_algebric_distance(vectors[i],vectors[j]) for i in range(len(vectors))] for j in range(len(vectors))]
    #is_collected:list = [False for _ in vectors]



def get_memory() -> float:
    return Process(getpid()).memory_info().rss / 1024 ** 2


def clean_fasta(input: str) -> list[str]:
    with open(input, 'r') as reader:
        return [l for l in reader if l[0] not in ['\n', '>']]


def compress_naive(input: str, output: str):
    init_memory = get_memory()
    lines = sorted(clean_fasta(input))
    max_memory = get_memory()
    with open(output, 'w') as writer:
        writer.write(''.join(lines))
    return abs(max_memory - init_memory)

# mettre la pitite fonction ici :)


if __name__ == "__main__":
    """
    parser = ArgumentParser()
    parser.add_argument("-i", "--input",
                        help="Path to a fasta-like file to re-order", type=str)
    parser.add_argument("-o", "--output",
                        help="Path to a txt file, output of the program", type=str)
    args = parser.parse_args()
    """
    kmers = [''.join(choice(['A', 'T', 'C', 'G'])
                     for _ in range(100)) for _ in range(1000000)]
    hashes = {k: hash_kmer(k) for k in kmers}
    print(len(hashes.keys()) == len(set(l for l in hashes.values())))
