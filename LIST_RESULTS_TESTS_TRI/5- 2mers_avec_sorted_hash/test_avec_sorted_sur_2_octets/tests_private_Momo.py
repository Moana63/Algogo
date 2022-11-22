from time import monotonic
from datetime import timedelta
from main_proba_2mers import compress_naive, write_file
from os import listdir, system, path
import matplotlib.pyplot as plt
from typing import Callable

FILES: list[str] = sorted([f for f in listdir(
    "Ecoli_100Kb") if f.split('.')[-1] == 'fasta'])
FILES_NO_EXT: list[str] = [f.split('.')[0] for f in FILES]
INPUT: list[str] = [f"Ecoli_100Kb/{p}" for p in FILES]
OUTPUT: list[str] = [f"{inp.split('.')[0]}.txt" for inp in INPUT]
FUNC_TEST: list[Callable] = [
    (compress_naive, {}), (write_file, {})]


def timer(func):
    """
    Decorator ; returns execution time of decorated func
    * arg : descrptor name of job
    """
    def wrapper(*args, **kwargs):
        start_time = monotonic()
        info = func(*args, **kwargs)
        end_time = monotonic()
        return {'time': timedelta(seconds=end_time - start_time).microseconds, **info}
    return wrapper


@timer
def my_tester(func: Callable, input: str, output: str, kwargs: dict) -> None:
    return {'memory': func(input, output, **kwargs), 'file': output}


if __name__ == "__main__":
    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, nrows=1, figsize=(15, 10))
    ax1.title.set_text('Fig. A : Sort function time spent')
    ax2.title.set_text('Fig. B : Sort function memory usage')
    ax3.title.set_text('Fig. C : Compression efficiency')

    for (func, kwargs) in FUNC_TEST:
        func_ret: list[tuple] = [
            my_tester(func, INPUT[i], OUTPUT[i], kwargs) for i, _ in enumerate(INPUT)]
        [system(f"gzip -f {out}") for out in OUTPUT]
        x: list = [ret['file'].split('.')[0].split('_')[-1]
                   for ret in func_ret]
        clean_files: list = [ret['file'].split('.')[0] for ret in func_ret]
        serie_time: list = [ret['time'] for ret in func_ret]
        serie_memory: list = [ret['memory'] for ret in func_ret]
        serie_du: list = [(path.getsize(
            f"{f}.txt.gz")/path.getsize(f"{f}.fasta.gz"))*100 for f in clean_files]
        ax1.plot(x, serie_time, label=f"{func.__name__}, {kwargs}")
        ax2.plot(x, serie_memory, label=f"{func.__name__}, {kwargs}")
        ax3.plot(x, serie_du, label=f"{func.__name__}, {kwargs}")

    handles, labels = ax1.get_legend_handles_labels()
    ax1.set_ylabel("Time (in microseconds)")
    ax2.set_ylabel("Peak memory used (Mb)")
    ax3.set_ylabel("Percentage of original file size")
    fig.legend(handles, labels, loc='upper center')
    plt.savefig("test.png")
    plt.show()
