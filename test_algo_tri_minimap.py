from itertools import groupby
# https://academic.oup.com/bioinformatics/article/34/18/3094/4994778?login=false#121258097


def hash(seq):
    hash_dict = {"A": 1, "T": 2, "C": 3, "G": 4, "N": 5}
    base = len(hash_dict)
    hash_val = sum([hash_dict[letter]*(base**i)
                   for i, letter in enumerate(seq)])
    return hash_val


def minimap(seq, len_window=5, kmer=6):
    list_kmer_window = [0 for _ in range(len_window)]
    minimiser_list = list()
    for i in range(len(seq)-kmer-1):
        for j in range(len_window):
            list_kmer_window[j] = (seq[i+j: i+j+kmer], i)
        min = sorted(list_kmer_window)[0]
        if min[0] not in minimiser_list:
            minimiser_list.append(min)
    return minimiser_list


# print(minimap("GCATGGGGCAAAAATGCCGAAGATGCGGTGCATAACGCCATCGTGCTGGAAGAGGTCGCTTATATGGGGATATTCTGCCGTCAGTTAGCGCCGCAGTTAC"))


def mini_index(ref_seq):
    dict_minimap = dict()
    minimiser_list = minimap(ref_seq.strip())
    for (minimiser, i) in minimiser_list:
        hash_minimiser = hash(minimiser)
        if hash_minimiser not in dict_minimap.keys():
            dict_minimap[hash_minimiser] = [i]
        else:
            dict_minimap[hash_minimiser] += [i]
    return dict_minimap


print(mini_index("GCATGGGGCAAAAATGCCGAAGATGCGGTGCATAACGCCATCGTGCTGGAAGAGGTCGCTTATATGGGGATATTCTGCCGTCAGTTAGCGCCGCAGTTAC"))

'''
def mini_index_multiple(list_ref_seq):
    dict_minimap = dict()
    for seq in list_ref_seq:
        minimiser_list = minimap(seq.strip())
        for (minimiser, i) in minimiser_list:
            hash_minimiser = hash(minimiser)
            if hash_minimiser not in dict_minimap.keys():
                dict_minimap[hash_minimiser] = [i]
            else:
                dict_minimap[hash_minimiser] += [i]
    return dict_minimap


get the minimiser's list (value) for each seq (key)
def mini_index(input):
    dict_minimap = dict()
    with open(input) as f:
        for index_line, line in enumerate(f):
            minimiser_list = minimap(line.strip())
            dict_minimap[index_line] = sorted([
                hash(minimiser) for minimiser in minimiser_list])
    return dict_minimap    



def ordering(input):
    dict_minimap = mini_index(input)
    score_list = 0
    line_closest_seq = 0
    dict_closest_lines = dict()
    for key in dict_minimap.keys():
        for index_line, value in enumerate(dict_minimap.values()):
            seq_ref = dict_minimap[key]
            new_score_list = intersection(seq_ref, value)
            if new_score_list > score_list and index_line != key:
                score_list = new_score_list
                line_closest_seq = index_line
        dict_closest_lines[key] = line_closest_seq
    return dict_closest_lines


def intersection(list_ref, list_2):
    # revoir calcul de score
    in_common = len([0 for value in list_ref if value in list_2])
    return in_common
'''
