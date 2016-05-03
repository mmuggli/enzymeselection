#!/usr/bin/python3

from itertools import *
from math import *

# simulate intervals from integer (i.e. discrete), semi-open ranges
two_intervals = [set(range(2,7 + 1))    | set(range(22, 27 + 1)),
                 set(range(0,4 + 1))    | set(range(25, 28 + 1)),
                 set(range(6, 10 + 1))  | set(range(19, 23 + 1)),
                 set(range(9, 13 + 1))  | set(range(16, 21 + 1)),
                 set(range(11, 18 + 1)) | set(range(29, 32 + 1))]

print(two_intervals)

vertices = list(range(len(two_intervals)))

all_pairs = [sorted((u, v)) for (u, v) in combinations(vertices, 2)]
#edges = [(u,v) for (u,v) in all_pairs if two_intervals[u] & two_intervals[v]]

    
#print(edges)

num_enzymes = len(all_pairs) + len(vertices)
num_bits = ceil( log2(num_enzymes) )
f_str = "{:0" + str(num_bits) + "b}"

enzymes = [f_str.format(e_num) for e_num in range(num_enzymes)]

def enz_str_of_node(v):
    return enzymes[v + len(all_pairs)]

def enz_str_of_pair(u, v):
    i = all_pairs.index(sorted((u,v)))
    print("index:", i,enzymes[i])
    return enzymes[i]

contigs = []

l = 0
for (u,v) in all_pairs:
    if two_intervals[u] & two_intervals[v]:
        contigs.append("#" * (l) + enz_str_of_pair(u,v))
        contigs.append("#" * (l + 1) + enz_str_of_pair(u,v))
        l += 1
    else:
        contigs.append("#" * l + enz_str_of_pair(u,v))
        contigs.append("#" * l + enz_str_of_pair(u,v))

print(contigs)        
