#!/usr/bin/python3

from itertools import *
from math import *

# simulate intervals from integer (i.e. discrete), semi-open ranges
# two_intervals = [set(range(2,7 + 1))    | set(range(22, 27 + 1)),
#                  set(range(0,4 + 1))    | set(range(25, 28 + 1)),
#                  set(range(6, 10 + 1))  | set(range(19, 23 + 1)),
#                  set(range(9, 13 + 1))  | set(range(16, 21 + 1)),
#                  set(range(11, 18 + 1)) | set(range(29, 32 + 1))]
#print('2-interval family from Fig. 1 of "On the parameterized complexity of multiple-interval graph problems" Michael R. Fellows a , Danny Hermelinb, Frances Rosamonda , Stephane Vialette c')

two_intervals = [set(range(0, 2 + 1)) | set(range(8,  10 + 1)),
                 set(range(1, 3 + 1)) | set(range(13, 15 + 1)),
                 set(range(4, 7 + 1)) | set(range(11, 14 + 1)),
                 set(range(5, 6 + 1)) | set(range(9,  12 + 1))]

#print("integer (discrete) transcription:", two_intervals)

vertices = list(range(len(two_intervals)))

all_pairs = [sorted((u, v)) for (u, v) in combinations(vertices, 2)]
#edges = [(u,v) for (u,v) in all_pairs if two_intervals[u] & two_intervals[v]]
antiedges = [(u,v) for (u,v) in all_pairs if not (two_intervals[u] & two_intervals[v])]
    
#print(edges)

num_enzymes = len(all_pairs) + len(vertices)
num_bits = ceil( log2(num_enzymes) )
f_str = "{:0" + str(num_bits) + "b}"

enzymes = [f_str.format(e_num) for e_num in range(num_enzymes)]
#print("edge enzymes in E:", enzymes[:len(all_pairs)])
#print("vertex enzymes in E:", enzymes[len(all_pairs):])
def enz_str_of_node(v):
    return enzymes[v + len(all_pairs)]

def enz_str_of_pair(u, v):
    i = all_pairs.index(sorted((u, v)))
    return enzymes[i]
# print("disjoint pairs:")
# for (u, v) in antiedges:
#     print(u, "(", enz_str_of_node(u), ") <-> ", v, "(",enz_str_of_node(v), ")  =", enz_str_of_pair(u, v))

econtigs = []


def dump(y, label, intervals):
    print(r"\node [right] at (0, " + str(y) + ") {" + label + "};")
    x = 3
    enzwid = .25
    step = .25
    for i in intervals:
        print("\draw [thick] (" +  str(x) + "," + str(y) + ") -- (" + str(x + i*step) + ", " + str(y) + ");")
        x += i * step
        print("\draw [thick, lightgray] (" +  str(x) + "," + str(y) + ") -- (" + str(x + enzwid) + ", " + str(y) + ");")
        x += enzwid

y = 100

ydec = .5
dump(y, "Edge selection:", [])
y -= ydec
# enzyme selection contigs
l = 0
for (u, v) in all_pairs:
    if two_intervals[u] & two_intervals[v]:
        dump(y, "("+str(u) +","+str(v)+")", [l])
        y-= ydec
        dump(y, "("+str(u) +","+str(v)+")", [l])
        y-=ydec
        #econtigs.append('#' * l + enz_str_of_pair(u, v))
        #econtigs.append('#' * l + enz_str_of_pair(u, v))
        l += 1
    else:
        # are these necessary? by construction, these will either not be in C' (if the edge enzyme is not selected) or else both be in C and by construction not contradict definition of C
        dump(y, "("+str(u) +","+str(v)+")", [l])
        y-= ydec
        dump(y, "("+str(u) +","+str(v)+")", [l+1])
        y-=ydec
        
        # econtigs.append('#' * (l) + enz_str_of_pair(u, v))
        # econtigs.append('#' * (l + 1) + enz_str_of_pair(u, v))
        l += 2
#print("\n".join(econtigs))
        
# (1) vertex selection: validateion of the selection of the endpoints of the edges
dump(y, "Vertex Selection - Endpoints of Edges:", [])
y -= ydec
v1acontigs = []
x = 1
lp = l #1 # l' el-prime
for (u, v) in antiedges:
    for k in range(len(antiedges) + 1):
        #v1acontigs.append('#' * x +  enz_str_of_pair(u, v) + '#' * lp + enz_str_of_node(u))
        dump(y, "("+str(u) +","+str(v)+"), "+str(u), [x, lp])
        y -= ydec
        lp += 1
    x += 1
#print("\n".join(v1acontigs))    

#dump(y, "Vertex selection contigs (1.b):", [])
#y -= ydec
v1bcontigs = []
x = 1
lpp = lp #1 # l' el-prime
for (u, v) in antiedges:
    for k in range(len(antiedges) + 1):
        v1bcontigs.append('#' * x +  enz_str_of_pair(u, v) + '#' * x +  enz_str_of_pair(u, v) + '#' * lpp + enz_str_of_node(v))
        dump(y, "("+str(u) +","+str(v)+"), ("+str(u) +","+str(v)+"), "+str(u), [x,x, lpp])        

        y -= ydec

        lpp += 1
    x += 1
#print("\n".join(v1bcontigs))
    
# (2) vertex selection: validateion of the selection of the edges of the endpoints
dump(y, "Vertex Selection - Edges of Endpoints:", [])
y -= ydec
v2acontigs = []
x = 1
lppp = lpp#1 # l' el-prime
for (u, v) in antiedges:
    for k in range(len(antiedges) + 1):
        #v2acontigs.append('#' * x +  enz_str_of_node(u) + '#' * lppp + enz_str_of_pair(u, v))
        dump(y,  str(u)+", ("+str(u) +","+str(v)+")", [x, lppp])
        y -= ydec
        lppp += 1
    x += 1
#print("\n".join(v2acontigs))

#dump(y, "Vertex selection contigs (2.b):", [])
#y -= ydec
v2bcontigs = []
x = 1
lpppp = lppp#1 # l' el-prime
for (u, v) in antiedges:
    for k in range(len(antiedges) + 1):
        v2bcontigs.append('#' * x +  enz_str_of_node(v) + '#' * lpppp + enz_str_of_pair(u, v))
        dump(y,  str(v)+", ("+str(u) +","+str(v)+")", [x, lpppp])
        y -= ydec
        
        lpppp += 1
    x += 1
#print("\n".join(v2bcontigs))

contigs = econtigs + v1acontigs + v1bcontigs + v2acontigs + v2bcontigs        
    
# code from http://stackoverflow.com/questions/3873361/finding-multiple-occurrences-of-a-string-within-a-string-in-python
# def findall(sub, string):
#     """
#     >>> text = "Allowed Hello Hollow"
#     >>> tuple(findall('ll', text))
#     (1, 10, 16)
#     """
#     index = 0 - len(sub)
#     try:
#         while True:
#             index = string.index(sub, index + len(sub))
#             yield index
#     except ValueError:
#         pass

# def digest(contig, enzymes):
#     positions = set()
#     for enzyme in enzymes:
#         positions |= set(findall(enzyme, contig))
#     return frozenset(positions)


# def valid(enzymes, contigs):
#     digested = [digest(contig, enzymes) for contig in contigs if digest(contig, enzymes)]
#     print("\n".join([str(s) for s in digested]))
#     return len(digested) == len(set(digested))

# selected_enzymes = ["0110", "1000", "0001"]
# print("valid enzymes(", selected_enzymes, "): ", valid(selected_enzymes, contigs))
# # for e in combinations(enzymes, 3):
# #     if valid(e, contigs): print(e)
