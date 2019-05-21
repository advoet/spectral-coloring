
# coding: utf-8

# In[1062]:


import graph_tool as gt
import graph_tool.generation as generation
import graph_tool.spectral as spectral
import graph_tool.topology as topology
import graph_tool.draw as draw


# In[1063]:


import matplotlib.pyplot as plt


# In[1064]:


import random
import itertools
import numpy
import numpy.linalg as linalg


# In[1065]:


# RANDOM GRAPH ON N VERTICES WITH DEGREE IN RANGE
def generateGraphByDegree(n, minDegree,  maxDegree):
    def deg():
        return random.randint(minDegree,maxDegree)

    g = generation.random_graph(n, deg, False)
    return g


# In[1127]:


draw.graph_draw(g)


# In[1126]:


g = generateGraphByProbability(20,.1)


# In[1116]:


for v in g.vertices():
    if v.out_degree() == 0:
        nonce = random.randint(1, g.num_vertices()-1)
        print(int(v), nonce)
        g.add_edge(int(v), int(v) + nonce )


# In[1128]:


# RANDOM GRAPH WITH EDGES GENERATED WITH FIXED PROBABILITY
def generateGraphByProbability(n, p):
    g = gt.Graph(directed=False)
    g.add_vertex(n)
    for (u,v) in itertools.combinations(g.get_vertices(), 2):
        if random.random() < p:
            g.add_edge(u,v)  
    # give the orphans a hand
    for v in g.vertices():
        if v.out_degree() == 0:
            nonce = random.randint(1, g.num_vertices()-1)
            g.add_edge(int(v), (int(v) + nonce) % g.num_vertices())
            
    return g


# In[1067]:


# ALEX

def rand_edge(g, u, v, p):
    """Adds an edge with probability p"""
    if random.random() < p:
        g.add_edge(g.vertex(u), g.vertex(v))

def generate_n_colorable(n, k, p):
    """Generates an n colorable graph with n*k nodes"""
    g = gt.Graph(directed=False)
    g.add_vertex(n*k)
    groups = [range(i*k,(i+1)*k) for i in range(0,n)]
    for vs in itertools.chain(*[itertools.product(range_1,range_2) for range_1,range_2 in itertools.combinations(groups, 2)]):
        rand_edge(g, vs[0], vs[1], p)
    return g

g = generate_n_colorable(20, 30, .4)


# In[1174]:


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


# In[1176]:


# ALEX 2

def generate_test_graph(n, k, a, c, m, bs):
    """Leighton's algorithm for generating graph with chromatic number"""
    g = gt.Graph(directed=False)
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
    x = random.randint(0,mod_m-1)
    for i in range(length):
        x = (a*x + c) % mod_m
        xs.append(x)
    return xs


# In[1165]:


g = generateGraphByProbability(50, .3)


# In[1167]:


g.get_vertices()


# In[1163]:


v = g.vertex(0)


# In[1162]:


v


# In[1083]:


int(v)


# In[1070]:


# STATS

def printStats(g):
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


# In[1144]:


# ADJACENCY SPECTRAL BOUNDS FUNCTION

def hoffmanQuotient(eigenvalues, m):
    num = sum(eigenvalues[:m])
    den = sum(eigenvalues[-m:])
    if den == 0:
        return 0
    return -num/den

def computeHoff(eigenvalues):
    return 1 + hoffmanQuotient(eigenvalues,1)
    
def computeGHoff(eigenvalues):
    return 1 + max([hoffmanQuotient(eigenvalues,m) for m in range(1,n)])
    
def computeWocElp(eigenvalues):
    posValues = list(filter(lambda x: x > 0, adjValues))
    negValues = list(filter(lambda x: x < 0, adjValues))
    Splus = sum(map(lambda x: pow(x,2), posValues))
    Sminus = sum(map(lambda x: pow(x,2), negValues))
    return 1 + Splus/Sminus

def computeWilf(eigenvalues):
    return 1 + eigenvalues[0]

def computeBounds(g):
    n = g.num_vertices()
    m = g.num_edges()
    average_degree = 2*m/n
    max_degree = max([vertex.out_degree() for vertex in g.vertices()])
    A = spectral.adjacency(g)
    Aw, Av = linalg.eig(A.todense())
    eigenvalues = sorted(Aw, reverse=True)
    
    results = []
    # compute eigen-bounds
    results.append(computeHoff(eigenvalues))
    results.append(computeGHoff(eigenvalues))
    results.append(computeWocElp(eigenvalues))
    results.append(computeWilf(eigenvalues))
    # stats
    results.append(max_degree)
    results.append(average_degree)
    # compute coloring
    results.append(getBestColoring(g)[0])
    
    return numpy.array(results)


# In[1145]:


def computeAverageBounds( generator, N):
    graphs = map(generator, range(N))
    bounds = map(computeBounds, graphs)
    averages = sum(bounds)/N
    result = zip(['hoff', 'ghoff', 'wocelp', 'wilf', 'max deg', 'avg deg','seq'], averages)
    return list(map(lambda x: x[0] + ': ' + str(x[1]), result))


# In[1157]:


# ALEX  3

def generate_strict_n_colorable(n, k, p):
    """Generates a graph with chromatic number n by introducing an n-clique"""
    g = generate_n_colorable(n, k, p)
    add_clique_2(n, g, [i*k for i in range(0,n)],0)
    
    return g

def generate_n_colorable(n, k, p):
    """Generates an n colorable graph with n*k nodes, edges with probability p"""
    g = gt.Graph(directed=False)
    g.add_vertex(n*k)

    #this line splits nodes into color groups, currently even split
    groups = [range(i*k,(i+1)*k) for i in range(0,n)]
    
    g.add_edge_list([pairs for pairs in itertools.chain(*[itertools.product(range_1,range_2) for range_1,range_2 in itertools.combinations(groups, 2)]) if random.random() < p])
    return g

def add_clique_2(clique_size, graph, ys, start_index):
    """Adds a clique to a graph from a vector ys of vertex indices"""
    graph.add_edge_list([pairs for pairs in itertools.combinations([ys[i % len(ys)] for i in range(start_index, start_index + clique_size)],2) if not graph.edge(pairs[0], pairs[1])])


# In[1170]:


print(generate_strict_n_colorable(10, 3, .5))


# In[1179]:


computeAverageBounds( lambda _: generate_test_graph(*tests[-1]), 5)


# In[1161]:


computeAverageBounds( lambda _: generate_strict_n_colorable(100, 3, .9), 5)


# In[1146]:


computeAverageBounds( lambda _: generateGraphByDegree(50, 2, 4), 15)


# In[1130]:


computeAverageBounds( lambda _: generateGraphByDegree(100, 2, 4), 15)


# In[1139]:


computeAverageBounds( lambda _: generateGraphByProbability(20, .5), 15)


# In[1135]:


computeAverageBounds( lambda _: generateGraphByProbability(100, .3), 15)


# In[1136]:


computeAverageBounds( lambda _: generateGraphByProbability(100, .5), 15)


# In[1137]:


computeAverageBounds( lambda _: generateGraphByProbability(100, .7), 15)


# In[1140]:


computeAverageBounds( lambda _: generateGraphByProbability(100, .9), 15)


# In[1143]:


computeAverageBounds( lambda _: generateGraphByProbability(20, .3), 15)


# In[1141]:


computeAverageBounds( lambda _: generateGraphByProbability(20, .1), 20)


# In[802]:


# NORMALIZED LAPLACIAN SPECTRAL BOUND(S)

B = spectral.laplacian(g, normalized=True)
Bw, Bv = linalg.eig(B.todense())

n = g.num_vertices()
m = g.num_edges()
average_degree = 2*m/n
lmax = max(Bw)
lBound = 1+1/(max(Bw)-1)
print('[?]: χ(g) ≥ ', lBound)

lBound = 1+max(Bw)/(max(Bw)-average_degree)
print('[Spielman(?)]: χ(g) ≥ ', lBound)


# In[959]:


def orderedColoring(g):
    prop = g.new_vertex_property("int32_t")
    randIndices = numpy.random.permutation(g.num_vertices())
    randDegrees = [g.vertex(i).out_degree() for i in randIndices]
    prop.a = [index for (degree, index) in sorted(zip(randDegrees, randIndices), key=lambda x: x[0],reverse=True)]
    degrees = [g.vertex(i).out_degree() for i in prop.a]
    colors = topology.sequential_vertex_coloring(g, order=prop)
    return (max(colors.a)+1,  colors)

def randomColoring(g):
    prop = g.new_vertex_property("int32_t")
    prop.a = numpy.random.permutation(g.num_vertices())
    colors = topology.sequential_vertex_coloring(g, order=prop)
    return (max(colors.a)+1,  colors)

def computeColoring(g, coloringType, iterations):
    colorings = list(map(lambda _: coloringType(g), range(iterations)))
    minIndex = numpy.argmin([coloring[0] for coloring in colorings])    
    
    return colorings[minIndex]

def getBestColoring(g):
    colorings = []
    colorings.append(computeColoring(g, randomColoring, 500))
    colorings.append(computeColoring(g, orderedColoring, 50))
    bestIndex = numpy.argmin([coloring[0] for coloring in  colorings])
    
    return colorings[bestIndex]


# In[957]:


c = getBestColoring(g)


# In[958]:


c


# In[941]:


# SKETCH WITH BEST SEQUENTIAL COLORING
draw.graph_draw(g, vertex_fill_color=c[1])


# In[1102]:


g = generateGraphByProbability(20, .1)


# In[1103]:


# PLAIN SKETCH

draw.graph_draw(g)


# In[ ]:


def colorings(k):
    for i in range(k):
        yield(i)


# In[ ]:


colorings(10)


# In[ ]:


import itertools


# In[ ]:


vMap = g.new_vertex_property("int32_t")
v = len(g.get_vertices())
iterator = itertools.product(range(3), repeat=v)


# In[ ]:


for thing in iterator:
    print(thing)


# In[ ]:


edges = g.edges()


# In[ ]:


# TERRIBLE  BRUTE FORCE ALGORITHM

colorable = False
index = 0
for coloring in iterator:
    vMap.a = list(coloring)
    success = True
    index += 1
    print(index)
    for e in g.edges():
        t = e.target()
        s = e.source()
        if vMap[t] == vMap[s]:
            success = False
            break
            
    if success:
        colorable = True
        break

if colorable:
    print('colorable')


# In[ ]:


coloring[e.target()]

