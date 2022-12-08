# FASTA-like file compression through reads reorganisation

## In short

Main purpose of this piece of software is to achieve better compression, by solely using data reorganisation, of a FASTA-like file. Across the process, we eventually lose all information apart from the lectures, as we do not keep headers in order to achieve better compression performance.

## Usage

In order to use the software, you have to provide a `.fasta` file and a output target at `.txt` format.
Then, you may specify a compression method. Default is `minimiser_presence_absence`, which achieves the best results across all tested scenarios.

```bash
read_organizer.py [-h] -i INPUT -o OUTPUT [--func [{minimisers_lexico,kmers_lexico,kmers_frequency,minimiser_presence_absence}]]
```



```bash
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to a fasta-like file to re-order
  -o OUTPUT, --output OUTPUT
                        Path to a txt file, output of the program
  --func [{minimisers_lexico,kmers_lexico,kmers_frequency,minimiser_presence_absence}]
                        Gives a method to sort reads.
```

Vérifier qu'on change pas la data

```bash
sort -u
diff
```