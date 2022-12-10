from sort_functions import kmers_frequency, kmers_lexico, minimiser_presence_absence, minimisers_lexico, clean_fasta
from time import process_time
from typing import Callable
from os import listdir, system, path
from numpy import asarray, mean, std
from memory_profiler import memory_usage
from argparse import ArgumentParser
import matplotlib.pyplot as plt

# reference files
FOLDER: str = "Ecoli_100Kb"
FILES: list[str] = sorted(
    [f for f in listdir(FOLDER) if f.split('.')[-1] == 'fasta'])
FILES_NO_EXT: list[str] = [f.split('.')[0] for f in FILES]
INPUT: list[str] = [f"{FOLDER}/{p}" for p in FILES]
OUTPUT: list[str] = [f"{inp.split('.')[0]}.txt" for inp in INPUT]
# list of tuples. each tuple contains => (func, {**kwargs})
FUNC_TEST: list[Callable] = [
    (minimisers_lexico, {'seed_size': i}) for i in range(2, 9)]
# function names for plot
FUNC_NAMES: list[str] = [f'ksize : {i}' for i in range(2, 9)]


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
def my_tester(func: Callable, input: str, output: str, kwargs: dict) -> dict:
    """aggregates results in a single dict

    Args:
        func (Callable): function to be executed
        input (str): file input
        output (str): file output
        kwargs (dict): all the keywords arguments for func

    Returns:
        dict: infos for plots
    """
    return {'memory': max(memory_usage((func, [input, output], kwargs))), 'file': output}


def compute_reference():
    """
    Creates a set of reference files for compression estimation
    """
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

    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, nrows=1, figsize=(40, 10))
    ax1.title.set_text('Fig. A : Sort function time spent')
    ax2.title.set_text('Fig. B : Sort function memory usage')
    ax3.title.set_text('Fig. C : Compression efficiency')

    for i, (func, kwargs) in enumerate(FUNC_TEST):  # iterating over all funcs
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
        # info if we work on a single reference file or multiple to change plot style
        single = False
        if len(asarray(series_time).transpose()) > 1:
            ax1.errorbar(x, [mean(serie) for serie in asarray(
                series_time).transpose()], yerr=[std(serie)/2 for serie in asarray(
                    series_time).transpose()], label=FUNC_NAMES[i], fmt='--o')
            ax2.errorbar(x,  [mean(serie) for serie in asarray(
                series_memory).transpose()], yerr=[std(serie)/2 for serie in asarray(
                    series_memory).transpose()], label=FUNC_NAMES[i], fmt='--o')
            ax3.errorbar(x,  [mean(serie) for serie in asarray(
                series_du).transpose()], yerr=[std(serie)/2 for serie in asarray(
                    series_du).transpose()], label=FUNC_NAMES[i], fmt='--o')
        else:
            single = True
            ax1.bar([FUNC_NAMES[i]], [mean(serie) for serie in asarray(series_time).transpose()], yerr=[std(
                serie)/2 for serie in asarray(series_time).transpose()], label=FUNC_NAMES[i], align='center', alpha=0.5, ecolor='black', capsize=10)
            ax2.bar([FUNC_NAMES[i]], [mean(serie) for serie in asarray(series_memory).transpose()], yerr=[std(
                serie)/2 for serie in asarray(series_memory).transpose()], label=FUNC_NAMES[i], align='center', alpha=0.5, ecolor='black', capsize=10)
            ax3.bar([FUNC_NAMES[i]], [mean(serie) for serie in asarray(series_du).transpose()], yerr=[std(
                serie)/2 for serie in asarray(series_du).transpose()], label=FUNC_NAMES[i], align='center', alpha=0.5, ecolor='black', capsize=10)

    ax1.set_ylabel("Time (in seconds)")
    ax2.set_ylabel("Peak memory used (Mb)")
    ax3.set_ylabel("Percentage of original file size")
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=len(labels))
    if single:
        ax1.tick_params(axis='x', which='both', bottom=False,
                        top=False, labelbottom=False)
        ax2.tick_params(axis='x', which='both', bottom=False,
                        top=False, labelbottom=False)
        ax3.tick_params(axis='x', which='both', bottom=False,
                        top=False, labelbottom=False)
    plt.savefig("report.png")
    plt.show()
