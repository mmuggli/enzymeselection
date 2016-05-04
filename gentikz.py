two_intervals = [((0,2) , (8,10)),
                 ((1,3) , (13,15)),
                 ((4,7) , (11,14)),
                 ((5,6) , (9,12))]

top = 5
left = 2

print(r"""
\documentclass{article}
\usepackage{tikz}
\begin{document}
\begin{tikzpicture}[scale=.5]
""")
for y, intervals in enumerate(two_intervals):
    print("\\node at (0, ", top - y, ") {$f_" + str( y),"$};")
    for start,end in intervals:
        print("\draw [thick] (", left + start * .5, ",", top - y, ") -- (", left + end * .5, ",", top - y, ");")
        print("\draw[fill] (",left + start * .5 ,",",top - y,") circle [radius=0.05];")
        print("\draw[fill] (",left + end * .5,",",top - y,") circle [radius=0.05];")

print(r"""        

\draw [line width=4,->] (12, 3.5) -- (14, 3.5);

\node at (18,5) {$v_0$};
\node at (21,5) {$v_1$};
\node at (21,2) {$v_2$};
\node at (18,2) {$v_3$};

\draw  (18,5) circle [radius=.5];
\draw  (21,5) circle [radius=.5];
\draw  (21,2) circle [radius=.5];
\draw  (18,2) circle [radius=.5];



\draw (18.5,5) -- (20.5,5);
\draw (21,4.5) -- (21,2.5);
\draw (20.5,2) -- (18.5,2);
\draw (18,2.5) -- (18,4.5);


\end{tikzpicture}  

\end{document}
""")
