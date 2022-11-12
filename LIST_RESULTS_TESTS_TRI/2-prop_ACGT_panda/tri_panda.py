# importing pandas package
import pandas as pd


def extract_proportions(filein):
    dict_prop = {}
    with open(filein) as f:
        list_seq = [seq.strip() for seq in f]
    for seq in list_seq:
        nb_A = seq.count("A")
        nb_C = seq.count("C")
        nb_G = seq.count("G")
        nb_T = seq.count("T")
        dict_prop[seq] = [nb_A, nb_C, nb_G, nb_T]
    df = pd.DataFrame.from_dict(
        dict_prop, orient='index', columns=["A", "C", "G", "T"])
    df.sort_values(["A", "C", "G", "T"], axis=0, ascending=[
                   False, False, False, False], inplace=True)
    return df
