# This file contains all the sorting methods that read_organizer uses in order to achieve its function
from collections import Counter
from itertools import product, chain

# list of all sort functions that can be accepted by the parser
PARSER_FUNCTIONS: list = [
    'minimisers_lexico',
    'kmers_lexico',
    'kmers_frequency',
    'minimiser_presence_absence'
]


def clean_fasta(input: str) -> list[str]:
    """Returns all the reads from a fasta-like file as a list

    Args:
        input (str): filepath to extract from

    Returns:
        list[str]: all the reads without their tags
    """
    with open(input, 'r') as reader:
        return [l.strip() for l in reader if l[0] not in ['\n', '>']]


def write_output(func) -> None:
    """Decorator to read input file as a list of reads, and writing out the returned list to the output

    Args:
        func (Callable): Targeted function to order reads. Must return a list
    """
    def wrapper(*args, **kwargs):
        with open(args[1], 'w') as writer:
            writer.write(
                '\n'.join(func(*args, **{'reads': clean_fasta(args[0]), **kwargs})))
        return None
    return wrapper

####################################################### INTERMEDIARY FUNCTIONS #######################################################


def frequency_kmer(read: str, seed_size) -> dict:
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


def indexation(list_seq: list, seed_size: int, len_read: int) -> dict:
    """Index all the read sequences from the fasta file according to their identifier.
    The identifier is a sequence of 0 and 1 linked to the proportions of différent kmers in the read sequence.

    Parameters
    ----------
    list_seq : list
        a list containing all the read sequence from the fasta file
    seed_size : int
        the size k of the kmer
    len_read : int
        the length of the read sequences in the file

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
        dico = frequency_kmer(read, seed_size)
        binary_seq = binary(dico, len_read, list_xmers, threshold)
        if binary_seq in index:
            index[binary_seq] += [i]
        else:
            index[binary_seq] = [i]
    return index


def frequency_minimizer(read: str, seed_size: int, len_window: int) -> dict:
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


def indexation_minimisers(list_seq: list, seed_size: int, len_window: int) -> dict:
    """Index all the read sequences from the fasta file according to their identifier. 
    The identifier is a sequence of 0 and 1 linked to the presence of different minimisers in the read sequence.

    Parameters
    ----------
    list_seq : list
        a list containing all the read sequence from the fasta file
    seed_size : int
       the size of the minimiser
    len_window : int
        size of the windo we seek minimizers in

    Returns
    -------
    dict
        an index with the identifier sequence of 0 and 1 as key and a list of the indexes in list_seq of the corresponding reads as values. 
    """
    index = dict()
    list_xmers = sorted(["".join(tuples)
                         for tuples in list(product('ATCG', repeat=seed_size))])
    for i, read in enumerate(list_seq):
        dico = frequency_minimizer(read, seed_size, len_window)
        binary_seq = binary_minimisers(dico, list_xmers)
        if binary_seq in index:
            index[binary_seq] += [i]
        else:
            index[binary_seq] = [i]
    return index

####################################################### SORT FUNCTIONS #######################################################


@write_output
def minimisers_lexico(input: str, output: str, reads: list = [], seed_size: int = 4, seed_number: int = 3, len_window: int = 33, **kwargs) -> list:
    """Sort a read file by the most present minimizers each read contains

    Args:
        input (str): input file
        output (str): output file
        seed_size (int, optional): size of minimizers. Defaults to 4.
        seed_number (int, optional): number of top common minimizers. Defaults to 3.
        len_window (int,optional): the length of the sliding window to extract the minimiser from the sequence, by default 33, which is the best option identified experimentally
    """
    return [reads[i] for i in [i for i, _ in sorted(enumerate([''.join(
        [key for key, _ in frequency_minimizer(read, seed_size, len_window).most_common(seed_number)]) for read in reads]), key=lambda x:x[1])]]


@write_output
def kmers_lexico(input: str, output: str, reads: list = [], seed_size: int = 4, seed_number: int = 3, **kwargs) -> list:
    """Sort a read file by the most present kmers each read contains

    Args:
        input (str): input file
        output (str): output file
        seed_size (int, optional): size of kmer. Defaults to 4.
        seed_number (int, optional): number of top common kmers. Defaults to 3.
    """
    return [reads[i] for i in [i for i, _ in sorted(enumerate([''.join(
        [key for key, _ in frequency_kmer(read, seed_size).most_common(seed_number)]) for read in reads]), key=lambda x:x[1])]]


@write_output
def kmers_frequency(input: str, output: str, reads: list = [], seed_size: int = 4, **kwargs) -> list:
    """Sort a read file by the kmers content of each read

    Parameters
    ----------
    input : str
        the file to sort
    output : str
        the sorted file
    reads : list, optional
        a list containing all the reads extracted from the input file, by default []
    seed_size : int, optional
        the size of kmer, used to sort the reads, by default 4, which is the best option identified experimentally

    Returns
    -------
    list
        the list containing the sorted reads

    """
    index = indexation(reads, seed_size, len(reads[0]))
    return list(chain(*[[reads[int(seq)] for seq in index[key]] for key in sorted(index.keys())]))


@write_output
def minimiser_presence_absence(input: str, output: str, reads: list = [], seed_size: int = 4, len_window: int = 33, **kwargs) -> list:
    """Sort a read file by the minimisers content of each read

    Parameters
    ----------
    input : str
        the file to sort
    output : str
        the sorted file
    reads : list, optional
        a list containing all the reads extracted from the input file, by default []
    seed_size : int, optional
        the size of kmer, used to sort the reads, by default 4, which is the best option identified experimentally
    len_window : int, optional
        the length of the sliding window to extract the minimiser from the sequence, by default 33, which is the best option identified experimentally

    Returns
    -------
    list
        the list containing the sorted reads
    """
    index = indexation_minimisers(reads, seed_size, len_window)
    return list(chain(*[[reads[int(seq)] for seq in index[key]] for key in sorted(index.keys())]))


if __name__ == "__main__":
    print("You shouldn't call this file directly, as it only contains functions. Try calling read_organiser.py instead.")
