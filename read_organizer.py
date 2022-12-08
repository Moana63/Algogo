import sort_functions
from argparse import ArgumentParser
from os import path, makedirs

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input",
                        help="Path to a fasta-like file to re-order",
                        type=str,
                        required=True)
    parser.add_argument("-o", "--output",
                        help="Path to a txt file, output of the program",
                        type=str,
                        required=True)
    parser.add_argument("--func",
                        default='minimiser_presence_absence',
                        const='minimiser_presence_absence',
                        nargs='?',
                        choices=['minimisers_lexico', 'kmers_lexico',
                                 'kmers_frequency', 'minimiser_presence_absence'],
                        help='Gives a method to sort reads.')
    args = parser.parse_args()

    # Verifiying if output is valid
    if args.output[-4:] != ".txt":
        print(TypeError(f"File {args.output} is not a valid .txt file."))
        exit()
    # Verifying if path to output is valid, creating it
    else:
        if '/' in args.output:
            makedirs(
                f"{'/'.join(args.output.split('/')[:-1])}/", exist_ok=True)

    if not path.exists(args.input) or not path.isfile(args.input):
        print(ValueError(f"Input file {args.input} does not exists."))
        exit()

    # Calling compress function by its name
    getattr(sort_functions, args.func)(args.input, args.output)
