import cProfile
from collections import Counter


a = Counter({'CATC': 3, 'CCGG': 3, 'CGGC': 3, 'ATCG': 2, 'TCGA': 2, 'ATCA': 2, 'CTTT': 2, 'TTTT': 2, 'ACCG': 2, 'GGCT': 2, 'GCTG': 2, 'CTGC': 2, 'GCTC': 2, 'CGAT': 1, 'GATA': 1, 'ATAA': 1, 'TAAA': 1, 'AAAC': 1, 'AACA': 1, 'ACAC': 1, 'CACA': 1, 'ACAT': 1, 'TCAA': 1, 'CAAA': 1, 'AAAA': 1, 'AAAT': 1, 'AATC': 1, 'TCAT': 1, 'ATCT': 1, 'TCTT': 1, 'TTTA': 1, 'TTAC': 1, 'TACC': 1, 'CCGC': 1, 'CGCA': 1, 'GCAT': 1, 'CGAG': 1, 'GAGG': 1, 'AGGC': 1, 'TGCT': 1, 'GCTT': 1, 'TTTG': 1, 'TTGC': 1, 'TGCA': 1, 'GCAC': 1, 'CACG': 1, 'ACGG': 1, 'CGGT': 1, 'GGTA': 1, 'GTAA': 1, 'TAAC': 1, 'AACG': 1, 'ACGC': 1, 'CGCC': 1, 'GCCT': 1, 'CCTG': 1, 'CTGT': 1, 'TGTT': 1, 'GTTT': 1, 'TTTC': 1, 'TTCC': 1, 'TCCC': 1, 'CCCG': 1, 'CTCC': 1, 'TCCG': 1, 'GGCC': 1, 'GCCA': 1, 'CCAG': 1, 'CAGC': 1, 'AGCT': 1, 'CTCA': 1, 'TCAC': 1, 'CACC': 1,
             'GGCG': 1, 'GCGT': 1, 'CGTC': 1, 'GTCG': 1, 'TCGC': 1, 'CGCT': 1, 'TGCC': 1, 'GCCC': 1})
b = Counter({'GATG': 3, 'ACGA': 3, 'ATGC': 2, 'TGCA': 2, 'GCAT': 2, 'AACG': 2, 'TGAT': 2, 'GATT': 2, 'ATTT': 2, 'GTAC': 2, 'CGAA': 2, 'GAAC': 2, 'TCTG': 2, 'CAAA': 2, 'ACTG': 2, 'CTGC': 2, 'TGCC': 2, 'CATA': 1, 'ATAA': 1, 'TAAC': 1, 'ACGC': 1, 'CGCG': 1, 'GCGG': 1, 'CGGA': 1, 'GGAT': 1, 'ATGT': 1, 'TGTG': 1, 'GTGA': 1, 'TTTT': 1, 'TTTC': 1, 'TTCG': 1, 'TCGC': 1, 'CGCC': 1, 'GCCG': 1, 'CCGT': 1, 'CGTC': 1, 'GTCG': 1, 'TCGG': 1, 'CGGG': 1, 'GGGG': 1, 'GGGT': 1, 'GGTA': 1, 'TACG': 1, 'CGAT': 1, 'TTTG': 1, 'TTGA': 1, 'ATGA': 1, 'TGAC': 1, 'GACC': 1, 'ACCG': 1, 'CCGA': 1, 'CGAC': 1, 'GACG': 1, 'AACA': 1, 'ACAA': 1, 'CAAT': 1, 'AATC': 1, 'ATCT': 1, 'CTGG': 1, 'TGGC': 1, 'GGCA': 1, 'GCAA': 1, 'AAAG': 1, 'AAGT': 1, 'AGTA': 1, 'TACT': 1, 'GCCC': 1, 'CCCA': 1, 'CCAA': 1, 'AAAT': 1, 'AATG': 1, 'GCCA': 1, 'CCAC': 1,
             'CACT': 1, 'CTGT': 1, 'TGTT': 1, 'GTTC': 1, 'TTCT': 1})

# kmer => groupes suivant les kmers les plus présents

cp = cProfile.Profile()

cp.enable()
for _ in range(100000):
    sum([(a[elt]-b[elt])**2 for elt in set().union(*a, *b)])
cp.print_stats()
cp.disable()

cp = cProfile.Profile()

cp.enable()
for _ in range(100000):
    keys = [*a.keys(), *b.keys()]
    sum([abs(a[elt]-b[elt]) for elt in keys])
cp.print_stats()
cp.disable()

cp = cProfile.Profile()

cp.enable()
for _ in range(100000):
    a.subtract(b)
    sum([abs(v) for v in a.values()])
cp.print_stats()
cp.disable()


cp = cProfile.Profile()

cp.enable()
for _ in range(100000):
    a.subtract(b)
    sum()
cp.print_stats()
cp.disable()
