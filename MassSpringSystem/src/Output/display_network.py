import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib.lines as mlines

#https://stackoverflow.com/questions/48242837/can-networkx-read-nodes-and-edges-from-different-files
#Answer by unutbu
edges=pd.read_csv('Results/bestMSEEdge.csv',sep=',')
G=nx.from_pandas_edgelist(edges,'from','to')

nodes = pd.read_csv('Results/bestMSENode.csv', sep=',')
data  = nodes.set_index('Node').to_dict('index').items()
G.add_nodes_from(data)

pos_x=nx.get_node_attributes(G,'x_pos')
pos_y=nx.get_node_attributes(G,'y_pos')

#https://stackoverflow.com/questions/5946236/how-to-merge-multiple-dicts-with-same-key-or-different-key?answertab=oldest#tab-top
#Answer by Eli Bendersky
pos = defaultdict(list)
for d in (pos_x, pos_y): # you can list as many input dicts as you want here
    for key, value in d.items():
        pos[key].append(value)

#https://stackoverflow.com/questions/62184321/set-node-color-based-on-attribute
#yatu
shape_fixed=[]
shape_other=[]
shape_dict= nx.get_node_attributes(G,'type')
for key, value in shape_dict.items():
    if(value=='FX'):
        shape_fixed.append(key)
    else:
        shape_other.append(key)

###############################################################################################################
#legend_map={'input node':'tab:green','feedback node':'plum','buckling':'yellow','fixed':'red','internal node':'white'}             
#f = plt.figure(1)
#ax = f.add_subplot(1,1,1)
#for label in legend_map:
#    ax.plot([0],[0],
#            color= legend_map[label],
#            label=label)
        
#https://stackoverflow.com/questions/47391702/matplotlib-making-a-colored-markers-legend-from-scratch
#tmdavison
input_node = mlines.Line2D([], [], color='tab:green', marker='o', linestyle='None',
                          markersize=10, label='input node')
feedback_node = mlines.Line2D([], [], color='plum', marker='o', linestyle='None',
                          markersize=10, label='feedback node')
buckling_node = mlines.Line2D([], [], color='yellow', marker='o', linestyle='None',
                          markersize=10, label='buckling node')
fixed_node = mlines.Line2D([], [], color='red', marker='s', linestyle='None',
                          markersize=10, label='fixed')
internal_node= mlines.Line2D([], [], color='lightgrey', marker='o', linestyle='None',
                          markersize=10, label='internal node')

################################################################################################################
#https://stackoverflow.com/questions/44801216/is-it-possible-to-mix-different-shaped-nodes-in-a-networkx-graph    
#https://groups.google.com/g/networkx-discuss/c/SinGFATZLaw
colour_map={'I':'tab:green','FB':'plum','B':'yellow','FX':'red','N':'lightgrey'}  
OG = G.subgraph(shape_other)
FG = G.subgraph(shape_fixed)
nx.draw_networkx_nodes(OG,pos,node_color=[colour_map[node[1]['type']]
                      for node in OG.nodes(data=True)],
                      node_shape='o',
                      edgecolors='black')
nx.draw_networkx_nodes(FG,pos,node_color=[colour_map[node[1]['type']]
                      for node in FG.nodes(data=True)],
                      node_shape='s',
                      edgecolors='black')
nx.draw_networkx_edges(G,pos)
plt.legend(handles=[input_node, feedback_node, buckling_node,fixed_node,internal_node]) 
plt.show() 


