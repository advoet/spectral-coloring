"""Script to generate graphs with given chromatic number"""
from graph_tool.all import *
from graph_tool.topology import is_bipartite
from random import randint, random
from itertools import product, chain, combinations
import numpy
import numpy.linalg as linalg
import time

tests = ((150, 5, 8401, 6859, 84035, (19,60,97,210)),
         (150, 5, 8401, 6859, 84035, (39,120,195,420)),
         (150,5,8401,6859,84035,(58,180,292,630)),
         (150,10,8401,6859,168070,(10,0,0,12,0,0,25,0,2,4)))

def main():
    g = Graph()
    g.add_vertex(200)
    start_1 = time.time()
    add_clique(100, g, range(0,200), 18)
    end_1 = time.time()
    print(end_1-start_1)

    g = Graph()
    g.add_vertex(200)
    start_1 = time.time()
    add_clique_2(100, g, range(0,200), 18)
    end_1 = time.time()
    print(end_1-start_1)

    g2 = generate_strict_three_colorable(5, .5)
    graph_draw(g2, vertex_text = g2.vertex_index, output="three.png", fmt="auto")
    g = generate_test_graph(*tests[3])
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

def generate_strict_n_colorable(n, k, p):
    """Generates a graph with chromatic number n by introducing an n-clique"""
    g = generate_n_colorable(n, k, p)
    add_clique_2(n, g, [i*k for i in range(0,n)])

def generate_n_colorable(n, k, p):
    """Generates an n colorable graph with n*k nodes, edges with probability p"""
    g = Graph(directed=False)
    g.add_vertex(n*k)

    #this line splits nodes into color groups, currently even split
    groups = [range(i*k,(i+1)*k) for i in range(0,n)]
    
    g.add_edge_list([pairs for pairs in chain(*[product(range_1,range_2) for range_1,range_2 in combinations(groups, 2)]) if random() < p])
    return g
    
def generate_test_graph(n, k, a, c, m, bs):
    """Leighton's algorithm for generating graph with chromatic number"""
    g = Graph(directed=False)
    g.add_vertex(n)
    
    xs = random_start_lcg(m, a, c, m)
    
    ys = [x % n for x in xs] #no longer uniform sequence
    
    start_index = 0
    clique_size = len(bs)+1
    for b in bs:
        for _ in range(0, b):
            add_clique_2(clique_size, g, ys, start_index)
            start_index = (start_index + clique_size) % n
        clique_size-=1
    return g

def random_start_lcg(length, a, c, mod_m):
    '''linear congruental generator'''
    xs = []
    x = randint(0,mod_m-1)
    for i in range(length):
        x = (a*x + c) % mod_m
        xs.append(x)
    return xs

def add_clique(clique_size, graph, ys, start_index):
    """Adds a clique to a graph from a vector ys of vertex indices

    replaced by add_clique_2
    """
    for i in range(0, clique_size-1):
        for j in range(i+1, clique_size):
            if graph.edge(graph.vertex(ys[(start_index + i) % len(ys)]),
                          graph.vertex(ys[(start_index + j) % len(ys)])):
                continue
            else:
                graph.add_edge(graph.vertex(ys[(start_index + i) % len(ys)]),
                           graph.vertex(ys[(start_index + j) % len(ys)]))

def add_clique_2(clique_size, graph, ys, start_index):
    """Adds a clique to a graph from a vector ys of vertex indices"""
    graph.add_edge_list([pairs for pairs in combinations([ys[i % len(ys)] for i in range(start_index, start_index + clique_size)],2) if not graph.edge(pairs[0], pairs[1])])
                
if __name__ == "__main__": main()
