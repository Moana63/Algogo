# This file contains all the sorting methods that read_organizer uses in order to achieve its function
from collections import Counter


def clean_fasta(input: str) -> list[str]:
    with open(input, 'r') as reader:
        return [l.strip() for l in reader if l[0] not in ['\n', '>']]


def write_output(output: str, reads_ordered: list) -> None:
    with open(output, 'w') as writer:
        writer.write('\n'.join(reads_ordered))


def frequency(read: str, seed_size) -> dict:
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


def compress_frequency_minimizer(input: str, output: str, ksize: int = 4, kmer_number: int = 3, len_window: int = 3) -> float:
    """Sort a read file by the most present kmers each read contains

    Args:
        input (str): input file
        output (str): output file
        ksize (int, optional): size of kmer. Defaults to 4.
        kmer_number (int, optional): number of top common kmers. Defaults to 3.

    Returns:
        float: _description_
    """
    reads: list = clean_fasta(input)
    ordered_reads: list = [reads[i] for i in [i for i, _ in sorted(enumerate([''.join(
        [key for key, _ in minimap_counter(read, len_window, ksize).most_common(kmer_number)]) for read in reads]), key=lambda x:x[1])]]
    write_output(output, ordered_reads)


def minimap_counter(seq, len_window, ksize):
    list_kmer_window = [0 for _ in range(len_window)]
    minimiser_list = list()
    for i in range(len(seq)-ksize-1):
        for j in range(len_window):
            list_kmer_window[j] = (seq[i+j: i+j+ksize], i)
        min = sorted(list_kmer_window)[0]
        if min[0] not in minimiser_list:
            minimiser_list.append(min)
    return Counter([minimizer for minimizer, idx in minimiser_list])


def compress_frequency_kmer(input: str, output: str, ksize: int = 4, kmer_number: int = 3) -> float:
    """Sort a read file by the most present kmers each read contains

    Args:
        input (str): input file
        output (str): output file
        ksize (int, optional): size of kmer. Defaults to 4.
        kmer_number (int, optional): number of top common kmers. Defaults to 3.

    Returns:
        float: _description_
    """
    reads: list = clean_fasta(input)
    ordered_reads: list = [reads[i] for i in [i for i, _ in sorted(enumerate([''.join(
        [key for key, _ in frequency(read, ksize).most_common(kmer_number)]) for read in reads]), key=lambda x:x[1])]]
    write_output(output, ordered_reads)


def minimisers_lexico(input: str, output: str) -> None:
    print("minim lexico")


def kmers_lexico(input: str, output: str) -> None:
    print("kmers lexico")


def kmers_frequency(input: str, output: str) -> None:
    print("kmers frequencies")


def minimiser_presence_absence(input: str, output: str) -> None:
    print("minimizers !!!!")


if __name__ == "__main__":
    print("You shouldn't call this file directly, as it only contains functions. Try calling read_organiser.py instead.")
