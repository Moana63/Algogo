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
                        choices=sort_functions.PARSER_FUNCTIONS,
                        help='Gives a method to sort reads.')
    parser.add_argument("--seed_size", default=4,
                        help='Defines a size for words we order by.')
    parser.add_argument("--len_window", default=33,
                        help='Defines a size for windows we go by.')
    parser.add_argument("--seed_number", default=3,
                        help='Defines a number of words we order by.')
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
    getattr(sort_functions, args.func)(args.input, args.output, kwargs={
        'seed_size': args.seed_size, 'len_window': args.len_window, 'seed_number': args.seed_number})
