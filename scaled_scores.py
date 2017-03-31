import networkx as nx
import sys

G = nx.Graph()

def get_max_incident(start, end):
   max=0
   try:
      for e in G.successors(start):
         if (e != end and G[start][e]['weight'] > max):
            max = G[start][e]['weight']
      for e in G.predecessors(start):
         if (e != end and G[e][start]['weight'] > max):
            max = G[start][e]['weight']
   except:
      for e in G.neighbors(start):
         if (e != end and G[start][e]['weight'] > max):
            max = G[start][e]['weight']

   return max

def get_max_weight(start,end):
   return max(get_max_incident(start, end), get_max_incident(end, start))

f = open('%s'% (sys.argv[1]))
for line in f:
   line = line.strip().split()

   #LNK tig00001320_pilon tig00019085_pilon E:B 103.138095238 21659 True True True
   #LNK tig00006867_pilon tig00009112_pilon E:E 0.214285714286 1 True True True
#   if line[0] != "LNK":
#     continue

#   edgeType=line[3].split(":")
   G.add_node(line[0])
   G.add_node(line[1])
   G.add_edge(line[0], line[1], weight=float(line[2]),isNeighbor="False", isGood="False")
f.close()

for (u, v, d) in G.edges(data='weight'):
   bestAlt=get_max_weight(u, v)
   if bestAlt==0:
      bestAlt=1

   print "%s %s %f %f %f %s %s"%(u, v, d, bestAlt, d/bestAlt, G[u][v]['isNeighbor'], G[u][v]['isGood'])
