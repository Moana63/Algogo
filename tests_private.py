from time import process_time
from datetime import timedelta
from main_Sisi import compress_naive, algebric_clustering, matrix_clustering, compress_by_kmer, clean_fasta
from typing import Callable
from os import listdir, system, path
import matplotlib.pyplot as plt
from numpy import asarray, transpose, mean, var, sqrt
from cProfile import Profile
from memory_profiler import memory_usage
from argparse import ArgumentParser

FILES: list[str] = sorted([f for f in listdir(
    "Ecoli_100Kb") if f.split('.')[-1] == 'fasta'])
FILES_NO_EXT: list[str] = [f.split('.')[0] for f in FILES]
INPUT: list[str] = [f"Ecoli_100Kb/{p}" for p in FILES]
OUTPUT: list[str] = [f"{inp.split('.')[0]}.txt" for inp in INPUT]
FUNC_TEST: list[Callable] = [(compress_by_kmer, {
                              'ksize': i, 'kmer_number': j}) for i in range(2, 6, 1) for j in range(3, 4, 1)]
# liste qui contient des tuples. chaque tuple contient => (le nom de la fonction, ses arguments dans un dictionnaire)


def timer(func):
    """
    Decorator ; returns execution time of decorated func
    * arg : descrptor name of job
    """
    def wrapper(*args, **kwargs):
        start_time = process_time()
        info = func(*args, **kwargs)
        end_time = process_time()
        return {'time': end_time - start_time, **info}
    return wrapper


@timer
def my_tester(func: Callable, input: str, output: str, kwargs: dict) -> None:
    return {'memory': max(memory_usage((func, [input, output], kwargs))), 'file': output}


def compute_reference():
    for inp in INPUT:
        reads: list = clean_fasta(inp)
        with open(f"reference/{inp.split('/')[-1]}", 'w') as writer:
            writer.write(''.join(reads))
    [system(f"gzip -f reference/{inp.split('/')[-1]}") for inp in INPUT]


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-p", "--passes",
                        help="Asks for a specific number of passes", type=int, default=10)
    parser.add_argument("-r", "--reference",
                        help="Tells if you need to create reference files", action='store_true')
    args = parser.parse_args()

    if args.reference:
        # to generate reference files at first execution
        system("mkdir reference/")
        compute_reference()

    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, nrows=1, figsize=(15, 10))
    ax1.title.set_text('Fig. A : Sort function time spent')
    ax2.title.set_text('Fig. B : Sort function memory usage')
    ax3.title.set_text('Fig. C : Compression efficiency')

    for (func, kwargs) in FUNC_TEST:
        series_time: list = []
        series_memory: list = []
        series_du: list = []
        for _ in range(int(args.passes)):  # averaging over 10 tries
            func_ret: list[dict] = [
                my_tester(func, INPUT[i], OUTPUT[i], kwargs) for i, _ in enumerate(INPUT)]
            [system(f"gzip -f {out}") for out in OUTPUT]
            x: list = [ret['file'].split('.')[0].split('_')[-1]
                       for ret in func_ret]
            clean_files: list = [ret['file'].split('.')[0] for ret in func_ret]
            series_time += [[ret['time'] for ret in func_ret]]
            series_memory += [[ret['memory'] for ret in func_ret]]
            series_du += [[(path.getsize(
                f"{f}.txt.gz")/path.getsize(f"reference/{f.split('/')[-1]}.fasta.gz"))*100 for f in clean_files]]
        ax1.errorbar(x, [mean(serie) for serie in asarray(
            series_time).transpose()], yerr=[sqrt(var(serie))/2 for serie in asarray(
                series_time).transpose()], label=f"{func.__name__}, {kwargs}", fmt='--o')
        ax2.errorbar(x,  [mean(serie) for serie in asarray(
            series_memory).transpose()], yerr=[sqrt(var(serie))/2 for serie in asarray(
                series_memory).transpose()], label=f"{func.__name__}, {kwargs}", fmt='--o')
        ax3.errorbar(x,  [mean(serie) for serie in asarray(
            series_du).transpose()], yerr=[sqrt(var(serie))/2 for serie in asarray(
                series_du).transpose()], label=f"{func.__name__}, {kwargs}", fmt='--o')

    handles, labels = ax1.get_legend_handles_labels()
    ax1.set_ylabel("Time (in seconds)")
    ax2.set_ylabel("Peak memory used (Mb)")
    ax3.set_ylabel("Percentage of original file size")
    fig.legend(handles, labels, loc='upper center')
    plt.savefig("test.png")
    plt.show()
