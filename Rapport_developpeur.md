Au format pdf
3 pages max

# Structure générale du code 
intro


# Fonctions générales

Communes à toutes les fonctions.

## Nettoyage des fichiers de séquençage : clean_fasta


## Ecriture des reads triés dans un fichier texte : @write_output

Le traitement des entrées/sorties s'effectue via un décorateur.

```python
def write_output(func) -> None:
    """Decorator to read input file as a list of reads, and writing out the returned list to the output

    Args:
        func (Callable): Targeted function to order reads. Must return a list
    """
```

On ajoute aux paramètres de `func` la liste des reads contenus dans le fichier (appel à `clean_fasta`)
L'écriture dans le fichier de sortie est appelé sur le retour de la fonction `func` décorée


# Première stratégie


## Algorithme de fréquence des kmers : kmers_lexico

### Fonction principale
```python
def kmers_lexico(input: str, output: str, reads: list = [], ksize: int = 4, kmer_number: int = 3) -> list:
    """Sort a read file by the most present kmers each read contains

    Args:
        input (str): input file
        output (str): output file
        ksize (int, optional): size of kmer. Defaults to 4.
        kmer_number (int, optional): number of top common kmers. Defaults to 3.
    """
```

### Fonctions annexes
```python

```

```python

```

```python

```


## Algorithme de fréquence des minimisers : minimisers_lexico

### Fonction principale
```python
def minimisers_lexico(input: str, output: str, reads: list = [], ksize: int = 4, kmer_number: int = 3, len_window: int = 3) -> list:
    """Sort a read file by the most present kmers each read contains

    Args:
        input (str): input file
        output (str): output file
        ksize (int, optional): size of kmer. Defaults to 4.
        kmer_number (int, optional): number of top common kmers. Defaults to 3.
    """
```

### Fonctions annexes
```python

```
```python

```
```python

```

# Seconde stratégie


## Algorithme de fréquence des kmers : kmers_frequency

### Fonction principale
```python

```

### Fonctions annexes
```python

```

```python

```

```python

```

## Algorithme de présence/absence des kmers : minimiser_presence_absence

### Fonction principale
```python

```

### Fonctions annexes
```python

```

```python

```

```python

```
