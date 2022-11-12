from itertools import groupby
# https://academic.oup.com/bioinformatics/article/34/18/3094/4994778?login=false#121258097


def hash(seq):
    hash_dict = {"A": 1, "T": 2, "C": 3, "G": 4, "N": 5}
    base = len(hash_dict)
    hash_val = sum([hash_dict[letter]*(base**i)
                   for i, letter in enumerate(seq)])
    return hash_val


def minimap(seq, len_window=8, kmer=3):
    list_kmer_window = [0 for _ in range(len_window)]
    minimiser_list = list()
    for i in range(len(seq)-kmer-len_window):
        for j in range(len_window):
            list_kmer_window[j] = seq[i+j: i+j+kmer]
        minimum = sorted(list_kmer_window)[0]
        if minimum not in minimiser_list:
            minimiser_list.append(minimum)
    return min(minimiser_list)


# renvoi un dictionnaire contenant
# comme clé = les minimisers minimaux de chaque séquence présente dans le fichier
# comme valeurs = la liste des séquences dans lesquelles ce minimiseur apparait


def mini_index(list_seq):
    dict_minimap = dict()
    for i, seq in enumerate(list_seq):
        minimiser = minimap(seq)
        hash_minimiser = hash(minimiser)
        if hash_minimiser not in dict_minimap.keys():
            dict_minimap[hash_minimiser] = [i]
        else:
            dict_minimap[hash_minimiser] += [i]
    return dict_minimap


def clustering_phase1(list_seq):
    dict_minimap = mini_index(list_seq)
    clusters = list()
    while dict_minimap:
        min_value = min((len(value), key)
                        for key, value in dict_minimap.items())
        list_values = dict_minimap[min_value[1]]
        del dict_minimap[min_value[1]]
        if list_values:
            clusters.append(list_values)
            for list_seq in dict_minimap.values():
                for seq in list_values:
                    if seq in list_seq:
                        list_seq.remove(seq)
    return clusters  # len(clusters), [len(value) for value in clusters]


def GC(ref_seq):
    return (ref_seq.count('G') + ref_seq.count('C')) / len(ref_seq)
