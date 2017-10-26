#Greedy agglomerative algorithm to create communities

import networkx as nx
import os
import operator as o
import math
import matplotlib.pyplot as plt


# calculating density of label 1
def a_1(g):
    sum = 0
    for i in nx.nodes_iter(g):
        length = len(g[i])
        if(g[i]['label']):
            sum = sum + (length - 1)**2
    return sum

#calculating density of label 0
def a_0(g):
    sum = 0
    for i in nx.nodes_iter(g):
        length = len(g[i])
        if (g[i]['label'] == 0):
            sum = sum + (length - 1) ** 2
    return sum

#calculating e for label 1
def e_1(g):
    count = 0
    for i in range(1,len(g)+1):
        for j in range(i+1, len(g)+1):
            if ((i,j) in nx.edges(g) and g[i]['label'] and g[j]['label']):
                count += 2
    return count

#calculating e for label 0
def e_0(g):
    count = 0
    for i in range(1,len(g)+1):
        for j in range(i+1, len(g)+1):
            if ((i,j) in nx.edges(g) and g[i]['label']==0 and g[j]['label']==0):
                count += 2
    return count

# calculating modularity for a graph with given partition with labels
def calcMod(g,l):
    e1 = 0.5*e_1(g)/l
    a1 = 0.25*a_1(g)/l**2
    sum_1 = e1 - a1
    sum_0 =  - ((a_0(g))/((2*l)**2))
    return sum_1 + sum_0

# labeling a node as 1
def addNode(g,i):
    if (g[i]['label'] == 0):
        g[i]['label']=1

#labeling node as 0
def removeNode(g,i):
    g[i]['label'] = 0

# inital matching for highest Q
def initModularity(g,totDegree):
    Q = calcMod(g,totDegree)
    qDiff = (-1)*math.inf
    qMap = {}

    for i in g.nodes_iter():
        addNode(g,i)
        maxQ = 0
        l = len(g[i]) - 1

        for j in (g[i]):
            if (l >= 1):
                addNode(g,j)
                qDiff1 = calcMod(g,totDegree) - Q
                if (qDiff1 > qDiff):
                    qDiff = qDiff1
                    maxQ = j
                removeNode(g,j)
                l = l - 1
        removeNode(g, i)
        if (maxQ):
            qMap[(i, maxQ)] = qDiff

    sorted(qMap.items(),key=o.itemgetter(1),reverse=True)
    (a,b) = qMap.pop()
    addNode(g,a)
    addNode(g,b)

    return g, (qDiff + Q), a, b

#Calculating Q as we continue the merge
def modularityForNode(g,totDeg,Q):
    qMap = {}
    qDiff = (-1)*math.inf
    maxQ = 0
    for i in g.nodes_iter():
        if (g[i]['label'] == 1):
            l = len(g[i]) - 1
            for j in g[i]:
                if (l >= 1 and g[j]['label']==0):
                    addNode(g,j)
                    qDiff1 = calcMod(g, totDeg) - Q
                    if (qDiff1 > qDiff):
                        qDiff = qDiff1
                        maxQ = j
                    removeNode(g, j)
                l = l - 1

            if (maxQ):
                qMap[(i,maxQ)] = qDiff
                if (g[maxQ]['label'] == 1):
                    removeNode(g,maxQ)


    qMap = sorted(qMap.items(),key=o.itemgetter(1),reverse=False)
    key, values = qMap.pop()
    c = key[0]
    d = key[1]

    addNode(g,d)
    addNode(g,c)
    return g, (Q+qDiff),qDiff, c , d

# calculating modularity for a graph
def totMod(g,totDeg):
    nodes = list(nx.nodes(g))
    qMergeMap = {}
    count = 1
    g, Q, a, b = initModularity(g, totDeg)
    addNode(g,a)
    addNode(g,b)
    #
    nodes.remove(a)
    nodes.remove(b)
    qMergeMap[1] = Q
    qDiff = 0

    while(nodes):
        count = count + 1
        g, Q, qDiff, a, b = modularityForNode(g,totDeg,Q)
        if (a in nodes):
            nodes.remove(a)
        if (b in nodes):
            nodes.remove(b)
        qMergeMap[count] = Q
        if(count == 17):
            print(nodes)
        print(count)
    return qMergeMap

def community_layout(g, partition):
    """
    Compute the layout for a modular graph.


    Arguments:
    ----------
    g -- networkx.Graph or networkx.DiGraph instance
        graph to plot

    partition -- dict mapping int node -> int community
        graph partitions


    Returns:
    --------
    pos -- dict mapping int node -> (float x, float y)
        node positions

    """

    pos_communities = _position_communities(g, partition, scale=3.)

    pos_nodes = _position_nodes(g, partition, scale=1.)

    # combine positions
    pos = dict()
    for node in g.nodes():
        pos[node] = pos_communities[node] + pos_nodes[node]

    return pos

def _position_communities(g, partition, **kwargs):

    # create a weighted graph, in which each node corresponds to a community,
    # and each edge weight to the number of edges between communities
    between_community_edges = _find_between_community_edges(g, partition)

    communities = set(partition.values())
    hypergraph = nx.DiGraph()
    hypergraph.add_nodes_from(communities)
    for (ci, cj), edges in between_community_edges.items():
        hypergraph.add_edge(ci, cj, weight=len(edges))

    # find layout for communities
    pos_communities = nx.spring_layout(hypergraph, **kwargs)

    # set node positions to position of community
    pos = dict()
    for node, community in partition.items():
        pos[node] = pos_communities[community]

    return pos

def _find_between_community_edges(g, partition):

    edges = dict()

    for (ni, nj) in g.edges():
        ci = partition[ni]
        cj = partition[nj]

        if ci != cj:
            try:
                edges[(ci, cj)] += [(ni, nj)]
            except KeyError:
                edges[(ci, cj)] = [(ni, nj)]

    return edges

def _position_nodes(g, partition, **kwargs):
    """
    Positions nodes within communities.
    """

    communities = dict()
    for node, community in partition.items():
        try:
            communities[community] += [node]
        except KeyError:
            communities[community] = [node]

    pos = dict()
    for ci, nodes in communities.items():
        subgraph = g.subgraph(nodes)
        pos_subgraph = nx.spring_layout(subgraph, **kwargs)
        pos.update(pos_subgraph)

    return pos


# reading data file
path = 'C:/Users/Sayali Sonawane/Documents/Fall-17/Network Analysis/PS3/zkcc-77/zkcc-77/karate_edges_77.txt'
fname = os.fsdecode(path)

#reading edgelist from file
g = nx.read_edgelist(fname, create_using=nx.DiGraph(), nodetype=int)
degree = nx.degree(g)
# number of edges
totDeg = sum(degree.values())/2
a = nx.adj_matrix(g)
for i in nx.nodes_iter(g):
    removeNode(g,i)

Q = totMod(g,totDeg)

#plotting Q vs number of merges
plt.xlabel('number of merges')
plt.ylabel('modularity Q')
xAxis = [0 for i in range(0,len(Q))]
yAxis = [0 for i in range(0,len(Q))]
l = len(xAxis)
for i in Q.keys():
    xAxis[i-1] = i
    yAxis[i-1] = Q[i]
plt.plot(xAxis,yAxis,label = 'Modularity')
plt.legend()
plt.show()

# plotting 2 communities of karate club
graph_pos = nx.spring_layout(g)
nx.draw_networkx_nodes(g,graph_pos,node_size=500,
                           alpha=0.3, node_color='blue')
nx.draw_networkx_edges(g,graph_pos,width=1,
                       alpha=0.3,edge_color='black')
nx.draw_networkx_labels(g, graph_pos,font_size=12,
                            font_family='sans-serif')
labels = {}
for j in nx.nodes(g):
    labels[j] = j
graph_pos = community_layout(g,nodeDict)
for com in set(nodeDict.values()):
    if (com == 0):
        color = 'blue'
    else:
        color = 'red'
    list_nodes = [nodes for nodes in nodeDict.keys() if nodeDict[nodes] == com]

    nx.draw_networkx_nodes(g, graph_pos, list_nodes, node_size=500, node_color=color, alpha=0.3)

nx.draw_networkx_edges(g,graph_pos, alpha=0.2, edge_color='blue')
nx.draw_networkx_labels(g,graph_pos,label=labels, font_size=16)

plt.show()
