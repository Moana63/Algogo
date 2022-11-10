from argparse import ArgumentParser
from os import getpid
from psutil import Process
from random import choice

k = 6

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

def indexing(input, k = k) -> None:
    index = dict()
    nb_lines = 0
    with open(input) as f:
        for line_index, line in enumerate(f):
            nb_lines += 1
            for i in range (0,len(line)-k,1):
                window = hash(line[i:i+k])
                if window in index:
                    index[window] += [(line_index, i)]
                else:
                    index[window] = [(line_index, i)]
    return nb_lines, index

def calc_score(input, k = k):
    line_index = 0 
    # get the index of the file
    len_file, index = indexing(input, k = k)
    print(len_file)
    # get the sequence of reference from the file (ref_seq)
    with open(input) as f:
        seq = "".join([ line for i, line in enumerate(f) if i == line_index ])
    # get the hash values of the kmers composing the sequence of reference (ref_seq)
    len_seq = len(seq)
    seq_hash = [(hash(seq[i:i+k]),i) for i in range (0,len(seq)-k,1)]
    # init score_list, 0 will be the score and i the index of the line associated with the score
    score_list = [[0,i] for i in range(len_file)]
    # for each kmer of the ref_seq, we get all the values associated in the index (kmer = key) and we calculate the score
    for hash_value in seq_hash:
        for value in index[hash_value[0]]:
    # the more the kmer is at a similar position between the ref_seq and the other sequence, the higher is the number added to the score
    # to obtain this we divide the length of the sequence by the difference between the index of the kmer on ref_seq and the index on the other sequence
    # Why => to avoid floatssssss
    # ATTENTION, NE PREND PAS EN COMPTE LE FAIT QUE LE MEME KMER PEUT ETRE REPETE PLUSIEURS FOIS DANS LA SEQUENCE DE REFERENCE => IMPACT SUR L'ALGO ?
            score_list[value[0]][0]+= len_seq//(abs(value[1]-hash_value[1])+1)
    # The list of scores is sorted by scores, from the highest to the lowest
    return sorted(score_list, reverse=True)

def write_ordered_kmers(input, output_file = "output.txt"):
    init_memory = get_memory()
    score_list = calc_score(input)
    with open(input) as f:
        seq_list = [line for line in f]
    max_memory = get_memory()
    with open(output_file, "w") as of:
        for item in score_list:
            of.write(seq_list[item[1]])
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
    write_ordered_kmers()
