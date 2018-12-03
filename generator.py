"""Script to generate graphs with given chromatic number"""
from graph_tool.all import *
from graph_tool.topology import is_bipartite
from random import randint, random
from itertools import product, chain, combinations

tests = ((5, 11, 150, 5, 8401, 6859, 84035, (19,60,97,210)),
         (5, 19, 150, 5, 8401, 6859, 84035, (39,120,195,420)),
         (5,24,150,5,8401,6859,84035,(58,180,292,630)),
         (10,11,150,10,8401,6859,168070,(10,0,0,12,0,0,25,0,2,4)))

def main():
    g2 = generate_strict_three_colorable(5, .5)
    graph_draw(g2, vertex_text = g2.vertex_index, output="three.png", fmt="auto")
    g = generate_test_graph(*tests[2])
    graph_draw(g, vertex_text = g.vertex_index, vertex_font_size=10,
               output_size = (5000,5000), output="test.png", fmt = "auto")

def generate_strict_three_colorable(k, p):
    """Generates a graph with 3k nodes of chromatic number 3"""
    g = generate_n_colorable(3, k, p)
    while is_bipartite(g):
        g = generate_n_colorable(3, k, p)
    return g

def generate_three_colorable(k, p):
    """Generates a three colorable graph with 3k nodes
    
    Not used anymore, replaced by generate_n_colorable with n=3
    """
    g = Graph(directed=False)
    g.add_vertex(3*k)
    for vs in chain(product(range(0,k),range(k,2*k)), product(range(0,k),range(2*k,3*k)), product(range(k,2*k),range(2*k,3*k))):
        rand_edge(g, vs[0], vs[1], p)
    return g

def rand_edge(g, u, v, p):
    """Adds an edge with probability p"""
    if random() < p:
        g.add_edge(g.vertex(u), g.vertex(v))

def generate_n_colorable(n, k, p):
    """Generates an n colorable graph with n*k nodes"""
    g = Graph(directed=False)
    g.add_vertex(n*k)
    groups = [range(i*k,(i+1)*k) for i in range(0,n)]
    for vs in chain(*[product(range_1,range_2) for range_1,range_2 in combinations(groups, 2)]):
        rand_edge(g, vs[0], vs[1], p)
    return g
    
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
    """Adds a clique to a graph from a vector ys of vertex indices"""
    for i in range(0, clique_size-1):
        for j in range(i+1, clique_size):
            if graph.edge(graph.vertex(ys[(start_index + i) % len(ys)]),
                          graph.vertex(ys[(start_index + j) % len(ys)])):
                continue
            else:
                graph.add_edge(graph.vertex(ys[(start_index + i) % len(ys)]),
                           graph.vertex(ys[(start_index + j) % len(ys)]))
                
if __name__ == "__main__": main()
