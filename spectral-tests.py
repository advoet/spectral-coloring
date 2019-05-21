
import graph_tool as gt
import graph_tool.generation as generation
import graph_tool.spectral as spectral
import graph_tool.topology as topology
import graph_tool.draw as draw
import matplotlib.pyplot as plt
import random
import itertools
import numpy
import numpy.linalg as linalg

from generator import *


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




#For a given graph g, compute the chromatic bounds from spectral theory and return
# list of results.
def compute_bounds(g):
    n = g.num_vertices()
    m = g.num_edges()
    average_degree = 2*m/n
    max_degree = max([vertex.out_degree() for vertex in g.vertices()])
    A = spectral.adjacency(g)
    Aw, Av = linalg.eig(A.todense())
    eigenvalues = sorted(Aw, reverse=True)
    
    results = []
    # compute eigen-bounds
    results.append(compute_hoff(eigenvalues))
    results.append(compute_gen_hoff(eigenvalues))
    results.append(compute_woc_elp(eigenvalues))
    results.append(compute_wilf(eigenvalues))
    # stats
    results.append(max_degree)
    results.append(average_degree)
    # compute coloring
    results.append(get_best_coloring(g)[0])
    
    return numpy.array(results)


#Uses a fixed graph generation, passed via lambda as generator, to generate N
#graphs and average the bounds from compute_bounds method
def compute_average_bounds( generator, N):
    graphs = map(generator, range(N))
    bounds = map(computeBounds, graphs)
    averages = sum(bounds)/N
    result = zip(['hoff', 'ghoff', 'wocelp', 'wilf', 'max deg', 'avg deg','seq'], averages)
    return list(map(lambda x: x[0] + ': ' + str(x[1]), result))

#Prints basic graph statistics
def print_stats(g):
    v = len(g.get_vertices())
    e = len(g.get_edges())
    average_degree = 2*e/v

    degrees = [vertex.out_degree() for vertex in g.vertices()]
    min_degree = min(degrees)
    max_degree = max(degrees)

    print('vertices: ', v)
    print('edges: ', e)
    print('min degree: ', min_degree)
    print('max degree: ', max_degree)
    print('average degree: ', average_degree)

############

# Functions to compute chromatic bounds based on the graph spectrum (passed
# as a list of eigenvalues)

############

def hoffman_quotient(eigenvalues, m):
    num = sum(eigenvalues[:m])
    den = sum(eigenvalues[-m:])
    if den == 0:
        return 0
    return -num/den

def compute_hoff(eigenvalues):
    return 1 + hoffman_quotient(eigenvalues,1)
    
def compute_gen_hoff(eigenvalues):
    return 1 + max([hoffman_quotient(eigenvalues,m) for m in range(1,n)])
    
def compute_woc_elp(eigenvalues):
    posValues = list(filter(lambda x: x > 0, adjValues))
    negValues = list(filter(lambda x: x < 0, adjValues))
    Splus = sum(map(lambda x: pow(x,2), posValues))
    Sminus = sum(map(lambda x: pow(x,2), negValues))
    return 1 + Splus/Sminus

def compute_wilf(eigenvalues):
    return 1 + eigenvalues[0]


###############
#
# A few coloring algorithms to test against spectral bounds
#
###############

def ordered_coloring(g):
    prop = g.new_vertex_property("int32_t")
    randIndices = numpy.random.permutation(g.num_vertices())
    randDegrees = [g.vertex(i).out_degree() for i in randIndices]
    prop.a = [index for (degree, index) in sorted(zip(randDegrees, randIndices), key=lambda x: x[0],reverse=True)]
    degrees = [g.vertex(i).out_degree() for i in prop.a]
    colors = topology.sequential_vertex_coloring(g, order=prop)
    return (max(colors.a)+1,  colors)

def random_coloring(g):
    prop = g.new_vertex_property("int32_t")
    prop.a = numpy.random.permutation(g.num_vertices())
    colors = topology.sequential_vertex_coloring(g, order=prop)
    return (max(colors.a)+1,  colors)

def compute_coloring(g, coloringType, iterations):
    colorings = list(map(lambda _: coloringType(g), range(iterations)))
    minIndex = numpy.argmin([coloring[0] for coloring in colorings])    
    
    return colorings[minIndex]

def get_best_coloring(g):
    colorings = []
    colorings.append(compute_coloring(g, random_coloring, 500))
    colorings.append(compute_coloring(g, ordered_coloring, 50))
    bestIndex = numpy.argmin([coloring[0] for coloring in  colorings])
    
    return colorings[bestIndex]



# compute_average_bounds( lambda _: generate_test_graph(*tests[-1]), 5)


# compute_average_bounds( lambda _: generate_strict_n_colorable(100, 3, .9), 5)


# compute_average_bounds( lambda _: generate_by_degree(50, 2, 4), 15)


# compute_average_bounds( lambda _: generate_by_degree(100, 2, 4), 15)


# compute_average_bounds( lambda _: generate_by_probability(20, .5), 15)


# compute_average_bounds( lambda _: generate_by_probability(100, .3), 15)


# compute_average_bounds( lambda _: generate_by_probability(100, .5), 15)


# compute_average_bounds( lambda _: generate_by_probability(100, .7), 15)


# compute_average_bounds( lambda _: generate_by_probability(100, .9), 15)


# compute_average_bounds( lambda _: generate_by_probability(20, .3), 15)


# compute_average_bounds( lambda _: generate_by_probability(20, .1), 20)



