"""Script to generate graphs with given chromatic number"""
from graph_tool.all import *
from random import randint

tests = ((5, 11, 150, 5, 8401, 6859, 84035, (19,60,97,210)),
         (5, 19, 150, 5, 8401, 6859, 84035, (39,120,195,420)),
         (5,24,150,5,8401,6859,84035,(58,180,292,630)),
         (10,11,150,10,8401,6859,168070,(10,0,0,12,0,0,25,0,2,4)))

def main():
    g = generate_test_graph(*tests[2])
    graph_draw(g, vertex_text = g.vertex_index, vertex_font_size=10,
               output_size = (5000,5000), output="test.png", fmt = "auto")

def generate_test_graph(chi, d, n, k, a, c, m, bs):
    g = Graph(directed=False)
    g.add_vertex(n)
    
    xs = [randint(0,m-1)]
    for i in range(0,n-1):
        xs.append((a*xs[i]+c) % m)
        
    ys = [x % n for x in xs]
    
    start_index = 0
    clique_size = len(bs)+1
    for i in range(0, len(bs)):
        for _ in range(0, bs[i]):
            add_clique(clique_size, g, ys, start_index)
            start_index = (start_index + clique_size) % n
        clique_size-=1
    return g

def add_clique(clique_size, graph, ys, start_index):
    for i in range(0, clique_size-1):
        for j in range(i+1, clique_size):
            if graph.edge(graph.vertex(ys[(start_index + i) % len(ys)]),
                          graph.vertex(ys[(start_index + j) % len(ys)])):
                continue
            else:
                graph.add_edge(graph.vertex(ys[(start_index + i) % len(ys)]),
                           graph.vertex(ys[(start_index + j) % len(ys)]))
                
if __name__ == "__main__": main()
