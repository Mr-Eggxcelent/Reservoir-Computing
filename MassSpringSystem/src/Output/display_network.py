import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib.lines as mlines
import matplotlib as mpl

mpl.rcParams['font.family'] = 'Avenir'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 2
#This script constructs the mass-spring network. Links to references I used are in the comments above each line of code.

#https://stackoverflow.com/questions/48242837/can-networkx-read-nodes-and-edges-from-different-files
#Answer by unutbu
edge_file="Results/quad_limit_cycles/bestMSEEdge.csv"
node_file="Results/quad_limit_cycles/bestMSENode_down.csv"

edges=pd.read_csv(edge_file,sep=',')
G=nx.from_pandas_edgelist(edges,'from','to')

nodes = pd.read_csv(node_file, sep=',')
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
input_node = mlines.Line2D([], [], color='#5fd35f', marker='o', linestyle='None',
                          markersize=10, label='input node')
feedback_node = mlines.Line2D([], [], color='#cd87de', marker='o', linestyle='None',
                          markersize=10, label='feedback node')
buckling_node = mlines.Line2D([], [], color='#5f8dd3', marker='o', linestyle='None',
                          markersize=10, label='buckling node')
fixed_node = mlines.Line2D([], [], color='red', marker='s', linestyle='None',
                          markersize=10, label='fixed')
internal_node= mlines.Line2D([], [], color='white', marker='o', linestyle='None',
                          markersize=10, label='internal node',markeredgecolor='black')

################################################################################################################
#https://stackoverflow.com/questions/44801216/is-it-possible-to-mix-different-shaped-nodes-in-a-networkx-graph    
#https://groups.google.com/g/networkx-discuss/c/SinGFATZLaw
colour_map={'I':'#5fd35f','FB':'#cd87de','B':'#5f8dd3','FX':'red','N':'white'}  
OG = G.subgraph(shape_other)
FG = G.subgraph(shape_fixed)
plt.figure(figsize=(6.4, 4.8)) 

nx.draw_networkx_nodes(OG,pos,node_color=[colour_map[node[1]['type']]
                      for node in OG.nodes(data=True)],
                      node_shape='o',
                      node_size=120,
                      edgecolors='black')
nx.draw_networkx_nodes(FG,pos,node_color=[colour_map[node[1]['type']]
                      for node in FG.nodes(data=True)],
                      node_shape='s',
                      node_size=120,
                      edgecolors='black')
nx.draw_networkx_edges(G,pos)
plt.legend(handles=[input_node,buckling_node,feedback_node,fixed_node,internal_node],loc="lower right") 
plt.show() 


