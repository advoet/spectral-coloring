"""Script to generate graphs with given chromatic number"""
from graph_tool.all import *
from graph_tool.topology import is_bipartite
import graph_tool.generation as generation
import graph_tool.spectral as spectral
import graph_tool.topology as topology
from random import randint, random
from itertools import product, chain, combinations
import numpy
import numpy.linalg as linalg
import time
import matplotlib.pyplot as plt


# Tests from Leighton '79' Appendix C
# (Number of nodes, chromatic number, a, c, m, (b))
tests = ((150, 5, 8401, 6859, 84035, (19,60,97,210)),
         (150, 5, 8401, 6859, 84035, (39,120,195,420)),
         (150, 5, 8401, 6859, 84035, (58,180,292,630)),
         (150, 10, 8401, 6859, 168070, (10,0,0,12,0,0,25,0,2,4)),
         (150, 10, 8401, 6859, 168070, (21, 0, 0, 24, 0, 0, 51, 0 ,48)),
         (150, 10, 8401, 6859, 168070, (31, 0, 0, 36, 0, 0, 76, 0, 72)),
         (150, 15, 8401, 6859, 252105, (4, 0, 0, 0, 0, 7, 0, 0, 0, 0, 12, 0, 0, 148)),
         (150, 15, 8401, 6859, 252105, (9, 0, 0, 0, 0, 15, 0, 0, 0, 0, 24, 0, 0, 297)),
         (150, 15, 8401, 6859, 252105, (13, 0, 0, 0, 0, 22, 0, 0, 0, 0, 36, 0, 0, 445)),
         (450, 5, 8401, 6859, 84034, (175, 540, 877, 1890)),
         (450, 5, 8401, 6859, 84034, (409, 1260, 2047, 4410)),
         (450, 15, 8401, 6859, 252105, (40, 0, 0, 0, 0, 67, 0, 0, 0, 0, 108, 0, 0, 1336)),
         (450, 15, 8401, 6859, 252105, (94, 0, 0, 0, 157, 0, 0, 0, 0, 252, 0, 0, 3118)),
         (450, 25, 8401, 6859, 420175, (13, 0, 0, 0, 0, 0, 0, 0, 0, 27, 0, 0, 0, 0, 0, 0, 0, 0, 40, 0, 0, 0, 0, 1336)),
         (450, 25, 8401, 6859, 420175, (31, 0, 0, 0, 0, 0, 0, 0, 0, 63, 0, 0, 0, 0, 0, 0, 0, 0, 94, 0, 0, 0, 0, 3118)))

def main():

    #Generate graph and image of a 3 colorable graph on 15=3*5 nodes
    g2 = generate_strict_three_colorable(5, .5)
    graph_draw(g2, vertex_text = g2.vertex_index, output="three.png", fmt="auto")

    #Unwrap test cases and generate image of test graph
    g = generate_test_graph(*tests[-1])
    graph_draw(g, vertex_text = g.vertex_index, vertex_font_size=10,
                output="test.png", fmt = "auto")


def generate_by_degree(n, minDegree,  maxDegree):
    """Generates a random graph on n vertices with degree in range(minDegree, maxDegree)

    Args:
        n (int): number of vertices
        minDegree (int): minimum degree required for each vertex
        maxDegree (int): maximum degree required for each vertex

    Returns:
        graph_tool graph
    """
    
    def deg():
        return random.randint(minDegree,maxDegree)

    g = generation.random_graph(n, deg, False)
    return g


def generate_by_edge_probability(n, p):
    """Generates a random graph on n vertices with edge probability p

    Args:
        n (int): number of vertices
        p : probability between 0 and 1 

    Returns: 
        graph_tool graph
    """

    g = gt.Graph(directed=False)
    g.add_vertex(n)
    for (u,v) in combinations(g.get_vertices(), 2):
        if random.random() < p:
            g.add_edge(u,v)  
    # give the orphans a hand
    for v in g.vertices():
        if v.out_degree() == 0:
            nonce = random.randint(1, g.num_vertices()-1)
            g.add_edge(int(v), (int(v) + nonce) % g.num_vertices())
            
    return g

def generate_strict_three_colorable(k, p):
    """Repeatedly applies method generate_n_colorable for n=3 until a non-bipartite graph is made.

    Args:
        k (int): number of vertices of each color, total number of vertices is 3*k
        p : probability between 0 and 1 of adding an edge between nodes of different colors

    Returns:
        graph_tool graph
    """

    g = generate_n_colorable(3, k, p)
    while is_bipartite(g):
        #If g happens to be 2 colorable, try again.
        g = generate_n_colorable(3, k, p)
    return g

def generate_three_colorable(k, p):
    """Generates a three colorable graph with 3k nodes, graph may be two-colorable (bipartite).

    Args:
        k (int): number of vertices of each color, total number of vertices is 3*k
        p : probability between 0 and 1 of adding an edge between nodes of different colors

    Returns:
        graph_tool graph
    """

    g = Graph(directed=False)
    g.add_vertex(3*k)
    for vs in chain(product(range(0,k),range(k,2*k)), product(range(0,k),range(2*k,3*k)), product(range(k,2*k),range(2*k,3*k))):
        rand_edge(g, vs[0], vs[1], p)
    return g

def rand_edge(g, u, v, p):
    """Introduces an edge in graph [g] between vertex [u] and vertex [v] with probability [p]

    Args:
        g (graph): graph tool graph, supporting g.vertex(u) and g.add_edge(vertex1, vertex2)
        u (int): vertex index
        v (int): vertex index
        p : probability between 0 and 1 of adding an edge
    """

    if random() < p:
        g.add_edge(g.vertex(u), g.vertex(v))

def generate_strict_n_colorable(n, k, p):
    """Generates a graph with chromatic number n by introducing an n-clique after running generate_n_colorable

    Args:
        n (int): chromatic number desired
        k (int): number of nodes of each color, there will be n*k vertices total
        p : probability between 0 and 1 of adding an edge between nodes of different colors

    Returns:
        graph_tool graph
    """

    #graph will be colorable with at most n colors
    g = generate_n_colorable(n, k, p)

    # graph is generated so that it can be colored with color groups [0,k), [k, 2*k), ..., [(n-1)*k, n*k)
    # We introduce a clique on nodes 0, k, 2k, ..., (n-1)*k so that any coloring requires at least n colors
    add_clique_2(n, g, [i*k for i in range(0,n)], 0)
    return g

def generate_n_colorable(n, k, p):
    """Generates an n colorable graph with n*k nodes, edges with probability p

    Args:
        n (int): output graph will be colorable with n colors
        k (int): number of nodes of each color, there will be n*k vertices total
        p : probability between 0 and 1 of adding an edge between nodes of different colors

    Returns:
        graph_tool graph
    """

    g = Graph(directed=False)
    g.add_vertex(n*k)

    #this line splits nodes into color groups, currently even split
    groups = [range(i*k,(i+1)*k) for i in range(0,n)]
    
    #adds an edge to the list between nodes from different groups (colors) with
    #probability p
    g.add_edge_list([pairs for pairs in chain(*[product(range_1,range_2) for range_1,range_2 in combinations(groups, 2)]) if random() < p])
    return g
    
def generate_test_graph(n, k, a, c, m, bs):
    """Uses Leighton's algorithm to generate a graph on [n] vertices with given chromatic number [k]

    Note the different input names from the rest of the methods to match variables from Leighton's paper. The purpose of this
    method is to create pseudo random graphs in a reproducable manner for repeated testing.

    Args:
        n (int): number of vertices
        k (int): chromatic number of output graph
        a (int): multiplier for lcg sequence, for any prime p and for p=4 we require p|m => p|(a-1)
        c (int): additive shift for lcg, positive and (c, m) = 1
        m (int): modulus for lcg, m>>n such that (n, m) = k and (c, m) = 1
        bs (int[]): integers (b_k, b_{k-1}, ... , b_2), where b_i is the number of cliques we add of size i. Require b_k > 0

    Returns:
        graph_tool graph
    """
    g = Graph(directed=False)
    g.add_vertex(n)
    
    xs = random_start_lcg(m, a, c, m)
    
    ys = [x % n for x in xs] #no longer uniform sequence
    
    #See Leighton '79
    start_index = 0
    clique_size = len(bs)+1
    for b in bs:
        for _ in range(0, b):
            add_clique_2(clique_size, g, ys, start_index)
            start_index = (start_index + clique_size) % n
        clique_size-=1
    return g

def random_start_lcg(length, a, c, mod_m):
    '''Linearly generates a sequence of numbers in Z/mZ

    Args:
        length (int): length of output sequence
        a (int): slope
        c (int): additive shift
        mod_m (int): modulus

    Returns:
        int[] array sequence of integers of given length
    '''
    
    xs = []
    x = randint(0,mod_m-1)
    for i in range(length):
        x = (a*x + c) % mod_m
        xs.append(x)
    return xs

def add_clique(clique_size, graph, ys, start_index):
    """Adds a clique to a graph from a vector ys of vertex indices

    replaced by add_clique_2, but this one is easier to read.
    """
    for i in range(0, clique_size-1):
        for j in range(i+1, clique_size):
            if graph.edge(graph.vertex(ys[(start_index + i) % len(ys)]),
                          graph.vertex(ys[(start_index + j) % len(ys)])):
                continue
            else:
                graph.add_edge(graph.vertex(ys[(start_index + i) % len(ys)]),
                           graph.vertex(ys[(start_index + j) % len(ys)]))

#Slightly optimized one line version of add_clique
def add_clique_2(clique_size, graph, ys, start_index):
    """Adds a clique to a graph from a vector ys of vertex indices

    introduces a clique on the vertices (ys[start_index], ys[start_index+1], ..., ys[start_index + clique_size - 1]) 
    wrapping around if necessary

    Args:
        clique_size (int): positive integer, size of clique to be added
        graph (graph_tool graph): its the graph
        ys (int[]): array of indices for distinct vertices in graph
        start_index (int): start position in array ys
    """

    graph.add_edge_list([pairs for pairs in combinations([ys[i % len(ys)] for i in range(start_index, start_index + clique_size)],2) if not graph.edge(pairs[0], pairs[1])])
                
if __name__ == "__main__": main()
