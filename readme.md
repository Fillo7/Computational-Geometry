Computational Geometry Project
==============================

Project is written in Processing: https://processing.org/download/

Controls:
A - Add new vertex under mouse cursor
D - Delete existing vertex under mouse cursor
Mouse click - Drag existing vertex to another place
R - Reload scene
P - Toggle vertex visibility
G - Generate 5 random vertices
C - Generate convex hull using Gift wrapping algorithm
V - Generate convex hull using Graham scan algorithm
T - Triangulation of polygon defined by all vertices in order of their addition, the last added vertex will be automatically
connected to the first, the polygon has to be y-monotone
T after C/V - Triangulation of generated convex hull
K - k-D tree, vertices located directly on the splitting edges belong to the upper left part
Y - Delaunay triangulation
Y after C/V - Delaunay triangulation of generated convex hull
U after Y - Voronoi diagram
