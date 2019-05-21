# spectral-coloring
CSE 599 Project: Spectral Graph Coloring

There are two files:

generator.py
contains methods for generating graphs. We provide methods for graph generation with several different requirements: edges exist with fixed probability, each node has degree within a fixed range, graph is n-colorable, graph has chromatic number n. We also implement Leighton's algorithm for reproducable pseudo-random test graphs based on a linear congruential generator.

spectral-tests.py
contains the code for testing various spectral bounds on the chromatic number of our graph. Featured are wilf's bound, a generalized version of the hoffman bound, and a conjectural bound due to Wocjan and Elphick.
