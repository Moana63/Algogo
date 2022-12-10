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
def kmers_frequency(input: str, output: str, reads: list = [], seed_size: int = 4) -> list:
    """Sort a read file by the kmers content of each read

    Parameters
    ----------
    input : str
        the file to sort
    output : str
        the sorted file
    reads : list, optional
        a list containing all the reads extracted from the input file, by default []
    seed_size : int, optional
        the size of kmer, used to sort the reads, by default 4, which is the best option identified experimentally

    Returns
    -------
    list
        the list containing the sorted reads

    """
```
Cette fonction récupère un dictionnaire indexant les reads à partir des fonctions annexes, l'index contient la métrique associé aux reads. Cette métrique est une séquence de 0 et de 1 rendant compte de façon simplifiée des proportions en kmers de la séquence. Les reads sont stockés via leur index dans la liste "reads" qui contient tous les reads du fichier, non triés. Cela permet d'économiser en mémoire. La fonction retourne ensuite une liste, contenant les reads triés par ordre alphanumérique. Elle est calculée de la façon suivante:

```python list(chain(*[[reads[int(seq)] for seq in index[key]] for key in sorted(index.keys())]))```

La boucle ```python for key in sorted(index.keys()``` tri les clés dans l'ordre alphanumérique. La boucle ```[[reads[int(seq)]for seq in index[key]]``` récupère les valeurs associées à chaque clé (leur position dans la liste reads) et renvoi les reads associés. On obtient à cette étape une liste contenant des listes. La méthode "chain*" permet d'applatir la liste en itérant avec "*" sur toutes les listes que contient la liste.
On pourrait choisir de trier les reads au moment de les écrire dans le fichier pour éviter de les stocker dans une liste, mais le choix a été fait d'avoir une seule fonction commune pour écrire le fichier trié et de faire moins d'action d'écriture dans le fichier d'output.

### Fonctions annexes
```python
def indexation(list_seq: list, seed_size: int, len_read: int = 100) -> dict:
    """Index all the read sequences from the fasta file according to their identifier.
    The identifier is a sequence of 0 and 1 linked to the proportions of différent kmers in the read sequence.

    Parameters
    ----------
    list_seq : list
        a list containing all the read sequence from the fasta file
    seed_size : int
        the size k of the kmer
    len_read : int, optional
        the length of the read sequences in the file, by default 100

    Returns
    -------
    dict
        an index with the identifier sequence of 0 and 1 as key and a list of the indexes in list_seq of the corresponding reads as values.
    """
```
Cette fonction génère le dictionnaire indexant les reads. Elle commence par récupérer un dictionnaire pour chaque read, contenant les kmers présents dans le read ainsi que leur nombre d'occurences. Elle simplifie ensuite le nombre d'occurence, en le passant à 1 si il est supèrieur à un seuil, ou à 0 sinon et récupérant le résultat sous la forme d'une séquence de 0 et de 1. Elle ajoute ensuite cette séquence dans le dictionnaire si elle n'existe pas encore et ajoute la position associée au read dans la liste des valeurs du dictionnaire.
On génère au début de la fonction "list_xmers" qui contient toutes les combinaisons de kmers possibles, de facon à s'en servir de référence pour que les positions de 0 et de 1 dans la séquence générée correspondent aux mêmes kmers pour tous les reads.

```python
def binary(dico: dict, len_read: int, list_xmers: list, threshold: float) -> str:
    """Take the proportions of each possible kmers in the sequence and returns a sequence of 0 and 1.
    For each kmer proportion, if the proportion is superior to the threshold, a 1 is added to the sequence. Else, a 0 is added.

    Parameters
    ----------
    dico : dict
        a dictionnary containing the proportions of each different kmers encountered in a sequence.
    len_read : int
        The length of the sequence
    list_xmers : list
        the list of all possible kmers of size k
    threshold : float
        the threshold above which the kmer proportion will be put to 1 and under which it will be put to 0.
    Returns
    -------
    str
        the sequence of 0 and 1 as a simplification of the proportions of kmers
    """
```
Cette fonction prend en entrée le dictionnaire contenant les occurences de kmers et permet de simplifier ce nombre d'occurences, en les passant à 1 si ils sont supèrieurs à un seuil, ou à 0 sinon. Elle renvoi le résultat sous la forme d'une séquence de 0 et de 1. On considère ici que si le nombre d'occurences du kmer est supèrieur au nombre d'occurence si les kmers étaient répartis da façon uniforme dans la séquence, alors ce kmer est surreprésenté et il sera associé à la valeure 1. Le calcul du seuil représente cette répartition uniforme théorique.

Utilise aussi la fonction frequency qui est commune aux 2 stratégies => CF

## Algorithme de présence/absence des kmers : minimiser_presence_absence

### Fonction principale
```python
def minimiser_presence_absence(input: str, output: str, reads: list = [], seed_size: int = 4, len_window: int = 33) -> list:
    """Sort a read file by the minimisers content of each read

    Parameters
    ----------
    input : str
        the file to sort
    output : str
        the sorted file
    reads : list, optional
        a list containing all the reads extracted from the input file, by default []
    seed_size : int, optional
        the size of kmer, used to sort the reads, by default 4, which is the best option identified experimentally
    len_window : int, optional
        the length of the sliding window to xtract the minimiser from the sequence, by default 33, which is the best option identified experimentally

    Returns
    -------
    list
        the list containing the sorted reads
    """
```
Cette fonction fonctionne de la même façon que la fonction principale 'kmers_frequency'. Le différence réside dans le fait que cette fois ci on récupère un nombre d'occurence de minimiseurs plutôt que de kmers et qu'on simplifie le problème en se basant sur la présence / absence de ces minimiseurs dans la séquence plutôt que leur fréquence d'apparition.  

### Fonctions annexes
```python
def indexation_minimisers(list_seq: list, seed_size: int, len_window: int = 33) -> dict:
    """Index all the read sequences from the fasta file according to their identifier. 
    The identifier is a sequence of 0 and 1 linked to the presence of different minimisers in the read sequence.

    Parameters
    ----------
    list_seq : list
        a list containing all the read sequence from the fasta file
    seed_size : int
       the size of the minimiser

    Returns
    -------
    dict
        an index with the identifier sequence of 0 and 1 as key and a list of the indexes in list_seq of the corresponding reads as values. 
    """
```
Cette fonction fonctionne sur le même principe que la fonction annexe 'indexation'.

```python
def binary_minimisers(dico: dict, list_xmers: list) -> str:
    """Take the proportions of each possible minimisers in the sequence and returns a sequence of 0 and 1. 
    For each minimisers proportion, if the minimiser is present, a 1 is added to the sequence. Else, a 0 is added.

    Parameters
    ----------
    dico : dict
        a dictionnary containing the proportions of each different minimisers encountered in a sequence.
    list_xmers : list
        the list of all possible minimisers

    Returns
    -------
    str
        the sequence of 0 and 1 as a simplification of the presence and absence of minimisers
    """
```
Cette fonction fonctionne sur le même principe que la fonction annexe 'binary', à la différence qu'on ne calcul pas de seuil, si le minimiseur est présent dans le read on ajoute un 1 à la séquence, sinon un 0.

```python
def frequency_minimizer(read: str, seed_size: int, len_window: int) -> dict:
    """Returns the frequency of minimisers per read

    Parameters
    ----------
    read : str
        a DNA read
    seed_size : int
        size of the minimiser
    len_window : int, optional
        length of the window sliding on the sequence, by default 10

    Returns
    -------
    dict
        a dictionnary containing the minimisers encountered in the sequence as key and their number of occurences as values
    """
```
Cette fonction fonctionne sur le même principe que la fonction annexe 'frequency', en incluant une fenêtre glissante qui parcourt la séquence et dont on peut ajuster la taille. 
Le minimiseur est récupéré en listant tous les kmers présents dans la fenêtre et en récupérant le plus petit (ordre lexicographique). On récupère aussi sa position j dans la séquence.
Pour améliorer la vitesse de parcourt de la séquence, le terme `python i += minimiser[1] + 1` afin d'accélérer le parcourt de la séquence par la fenêtre glissante. Dés que l'on a trouvé un minimiseur, on déplace la fenêtre de façon à dépasser ce minimiseur avant de recommencer à en chercher un. On utilise pour cela la position 'j' du minimiseur dans la séquence stocké dans `python minimiseur[1]`.
