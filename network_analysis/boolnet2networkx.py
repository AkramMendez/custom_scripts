# Este es un programa para generar redes bipartitas a partir de datos
# This script converts  network from boolnet format into a tabular #targets, factors format
# Returns some general graph measures using networkx

import math, sys
import networkx as nx
import matplotlib.pyplot as plt

# cargar archivo desde linea de comando
print sys.argv[1]
f = open(sys.argv[1], 'r')
G=nx.DiGraph()

# convertir a red
for line in f:
    if line != 'targets, factors\n':
    # quita caracteres inecesarios
        for char in line:
            if char in '&!|()\n':
                line = line.replace(char,'')
        # generar array
        temp = line.split(', ')
        temp[1] = temp[1].split()
        # declarar nodos
        G.add_node(temp[0]) #declare nodes
        for t in temp[1]: #declare edges
            G.add_edge(t, temp[0])
f.close()


#medidas estadisticas
#agregen las que faltan
print 'degree_centrality', nx.degree_centrality(G)
print 'eccentricity_undirected', nx.eccentricity(G.to_undirected())
print 'clustering_undirected', nx.clustering(G.to_undirected())


#Convert to other formats and draw
nx.write_gml(G, sys.argv[1]+'.gml') #exportar
nx.draw(G) #visualisar
plt.show()
