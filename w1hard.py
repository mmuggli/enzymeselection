#!/usr/bin/python3

from itertools import *
from math import *

# simulate intervals from integer (i.e. discrete), semi-open ranges
# two_intervals = [set(range(2,7 + 1))    | set(range(22, 27 + 1)),
#                  set(range(0,4 + 1))    | set(range(25, 28 + 1)),
#                  set(range(6, 10 + 1))  | set(range(19, 23 + 1)),
#                  set(range(9, 13 + 1))  | set(range(16, 21 + 1)),
#                  set(range(11, 18 + 1)) | set(range(29, 32 + 1))]
print('2-interval family from Fig. 1 of "On the parameterized complexity of multiple-interval graph problems" Michael R. Fellows a , Danny Hermelinb, Frances Rosamonda , Stephane Vialette c')

two_intervals = [set(range(0,2 + 1)) | set(range(8,10 + 1)),
                 set(range(1,3 + 1)) | set(range(13,15 + 1)),
                 set(range(4,7 + 1)) | set(range(11,14 + 1)),
                 set(range(5,6 + 1)) | set(range(9,12 + 1))]

print("integer (discrete) transcription:", two_intervals)

vertices = list(range(len(two_intervals)))

all_pairs = [sorted((u, v)) for (u, v) in combinations(vertices, 2)]
#edges = [(u,v) for (u,v) in all_pairs if two_intervals[u] & two_intervals[v]]
antiedges = [(u,v) for (u,v) in all_pairs if not (two_intervals[u] & two_intervals[v])]
    
#print(edges)

num_enzymes = len(all_pairs) + len(vertices)
num_bits = ceil( log2(num_enzymes) )
f_str = "{:0" + str(num_bits) + "b}"

enzymes = [f_str.format(e_num) for e_num in range(num_enzymes)]
print("edge enzymes in E:", enzymes[:len(all_pairs)])
print("vertex enzymes in E:", enzymes[len(all_pairs):])
def enz_str_of_node(v):
    return enzymes[v + len(all_pairs)]

def enz_str_of_pair(u, v):
    i = all_pairs.index(sorted((u, v)))
    return enzymes[i]
print("disjoint pairs:")
for (u, v) in antiedges:
    print(u, "(",enz_str_of_node(u),") <-> ",v, "(",enz_str_of_node(v), ")  =",enz_str_of_pair(u, v))

contigs = []
print("Enzyme selection contigs:")
# enzyme selection contigs
l = 0
for (u, v) in all_pairs:
    if two_intervals[u] & two_intervals[v]:
        print('#' * l + enz_str_of_pair(u, v))
        print('#' * l + enz_str_of_pair(u, v))
        l += 1
    else:
        # are these necessary? by construction, these will either not be in C' (if the edge enzyme is not selected) or else both be in C and by construction not contradict definition of C
        print('#' * (l) + enz_str_of_pair(u, v))
        print('#' * (l + 1) + enz_str_of_pair(u, v))
        l += 2

        
# (1) vertex selection: validateion of the selection of the endpoints of the edges
print("Vertex selection contigs (1.a):")
x = 1
lp = 1 # l' el-prime
for (u, v) in antiedges:
    for k in range(len(antiedges) + 1):
        print('#' * x +  enz_str_of_pair(u, v) + '#' * lp + enz_str_of_node(u))
        lp += 1
    x += 1

print("Vertex selection contigs (1.b):")
x = 1
lpp = 1 # l' el-prime
for (u, v) in antiedges:
    for k in range(len(antiedges) + 1):
        print('#' * x +  enz_str_of_pair(u, v) + '#' * x +  enz_str_of_pair(u, v) + '#' * lpp + enz_str_of_node(v))       
        lpp += 1
    x += 1

# (2) vertex selection: validateion of the selection of the edges of the endpoints
print("Vertex selection contigs (2.a):")
x = 1
lppp = 1 # l' el-prime
for (u, v) in antiedges:
    for k in range(len(antiedges) + 1):
        print('#' * x +  enz_str_of_node(u) + '#' * lppp + enz_str_of_pair(u, v))       
        lppp += 1
    x += 1


print("Vertex selection contigs (2.b):")
x = 1
lpppp = 1 # l' el-prime
for (u, v) in antiedges:
    for k in range(len(antiedges) + 1):
        print('#' * x +  enz_str_of_node(v) + '#' * lpppp + enz_str_of_pair(u, v))       
        lpppp += 1
    x += 1


        
    
print("\n".join(contigs))

