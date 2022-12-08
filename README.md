# FASTA-like file compression through reads reorganisation

## In short

Main purpose of this piece of software is to achieve better compression, by solely using data reorganisation, of a FASTA-like file. Across the process, we eventually lose all information apart from the lectures, as we do not keep headers in order to achieve better compression performance.

## Usage

In order to use the software, you have to provide a `.fasta` file and a output target at `.txt` format.
Then, you may specify a compression method. Default is `minimiser_presence_absence`, which achieves the best results across all tested scenarios.

```bash
read_organizer.py [-h] -i INPUT -o OUTPUT [--func [{minimisers_lexico,kmers_lexico,kmers_frequency,minimiser_presence_absence}]]
```

## Parameters

Command-line tool allows you to interact with the software with those arguments :

```bash
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to a fasta-like file to re-order
  -o OUTPUT, --output OUTPUT
                        Path to a txt file, output of the program
  --func [{minimisers_lexico,kmers_lexico,kmers_frequency,minimiser_presence_absence}]
                        Gives a method to sort reads. Default is 'minimiser_presence_absence'.
```

The `--func` argument allows you to pick a method to reorder the FASTA-like file.

### --func minimiser_presence_absence (default)

This method uses the presence or absence of minimizers, given a read-specific threshold, to reorder the file.

### --func minimisers_lexico

This method uses the lexicographical order of most common minimizers to reorder the file.

### --func kmers_frequency

This method uses the frequency of kmers, given a read-specific threshold, to reorder the file.

### --func kmers_lexico

This method uses the lexicographical order of most common kmers to reorder the file.