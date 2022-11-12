from argparse import ArgumentParser
from os import getpid
from psutil import Process
from random import choice
from test_algo_tri_minimap import clustering_phase1, GC


def hash_kmer(kmer: str) -> int:
    return sum(({'A': 1, 'C': 2, 'G': 3, 'T': 4, }.get(l, 0)*(5**i) for i, l in enumerate(kmer)))


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


def write_file(filein, fileout):
    init_memory = get_memory()
    list_seq = sorted(clean_fasta(filein))
    clusters = clustering_phase1(list_seq)
    max_memory = get_memory()
    with open(fileout, "w") as out_file:
        for group in clusters:
            group = sorted(group, key=lambda x: GC(list_seq[x]))
            for seq in group:
                out_file.write(list_seq[seq]+'\n')
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
    kmers = [''.join(choice(['A', 'T', 'C', 'G'])
                     for _ in range(100)) for _ in range(1000000)]
    hashes = {k: hash_kmer(k) for k in kmers}
    print(len(hashes.keys()) == len(set(l for l in hashes.values())))
    write_file()
