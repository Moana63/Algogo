Au format pdf
3 pages max

# Structure générale du code 
intro


# Fonctions générales

La fonction `frequency_minimizer` consiste à retourner, pour le read donné en argument, le dictionnaire contenant les minimiseurs de ksize de long d'une fenêtre glissante sur le read de len_window de long. Pour chaque position, on récupère un minimiseur, qui se trouve inclus au dictionnaire comptabilisant ceux-ci.
L'objet retourné est un dictionnaire, au format `minimiseur : nombre d'occurences du minimiseur` pour chaque minimiseur rencontré dans le read.

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
[
    reads[i] for i in
    [
        i for i, _ in sorted(
            enumerate(
                [
                    ''.join(
                        [
                            key for key, _ in frequency(read, ksize).most_common(kmer_number)
                        ]
                    ) for read in reads
                ]
            ), key=lambda x:x[1]
        )
    ]
]
```

On utilise ici une compréhension de liste qui, pour chaque read, extrait les kmer_number ksize-mers les plus communs.
Ensuite, on concatène grâce à la fonction join les kmers en une signature, que l'on trie lexicographiquement en fonction de cette signature. Enfin, par compréhension de liste, on trie les reads en fonction de l'index de la signature.

## Algorithme de fréquence des minimisers : minimisers_lexico

### Fonction principale

La seule différence par rapport à `kmers_lexico` est dans l'appel à la fonction permettant d'obtenir le comptage des minimiseurs au lieu des kmers. On fait ici appel à la fonction globale `frequency_minimizer` et non `frequency`.

```python
frequency_minimizer(read, ksize, len_window).most_common(kmer_number)
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
