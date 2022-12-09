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
    Args: func (Callable): Targeted function to order reads. Must return a list
    """
```

On ajoute aux paramètres de `func` la liste des reads contenus dans le fichier (appel à `clean_fasta`)
L'écriture dans le fichier de sortie est appelé sur le retour de la fonction `func` décorée.
Cela permet d'implémenter et de maintenir facilement les vérifications sur les fichiers d'entrée et de sortie, sans impacter la capacité à mesurer la mémoire et le temps requis pour l'exécution de l'algorithme de tri.

Cela permet d'implémenter toute autre nouvelle méthode de tri sans avoir à se questionner sur les entrées/sorties, tant que la fonction nouvellement implémentée prend en paramètres un fichier d'entrée, un de sortie, et un dictionnaire optionnel de reads. Toute nouvelle méthode doit suivre cette signature minimale :

```python
@write_output
def some_sort_function(input: str, output: str, reads: list = []) -> list:
    """
    Parameters
    ----------
    input : (str) the FASTA-like file to sort
    output : (str) the sorted file
    reads : (list, optional) reads extracted from the input file by the decorator.

    Returns
    -------
    (list) the list containing the sorted reads
    """
```


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
