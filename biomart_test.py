# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=2>

# Set up environment

# <codecell>

# general
import os
import sys
import urllib.request
from urllib.request import Request, urlopen
from urllib.parse import urlencode
import wget
import ssl # required for request to biomart on windows

# <codecell>

# analysis
import pandas as pd
import networkx as nx
import numpy as np
import scipy as sp
import scipy.spatial as spspat
import scipy.cluster as spclu

# <codecell>

# plotting
import matplotlib as mpl
import matplotlib.pyplot as plt

# <codecell>

# global vars
DIR = os.getcwd()
OUTPUT_DIR = DIR + '\\' + 'output\\'
if not os.path.exists(OUTPUT_DIR):
    print('Making ' + OUTPUT_DIR)
    os.makedirs(OUTPUT_DIR)
# PROXY = os.environ.get('http_proxy', 'http://eu-chbs-PROXY.eu.novartis.net:2011/')

# how long should Nmers be when parsing UTR sequences
NMER_LEN = 8

# <headingcell level=2>

# Biomart

# <headingcell level=5>

# Set up query

# <codecell>

# obtained from http://www.ensembl.org/biomart/
# enter parameters and retieve XML rather than results
# data request
query = """
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "ensembl_transcript_id" />
		<Attribute name = "3utr" />
	</Dataset>
</Query>
"""
# base url for request
baseurl = 'http://www.biomart.org/biomart/martservice'
# alt: http://www.ensembl.org/biomart/martview/d39bd0cb5716fe00e31e6fbd8a12f030?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.sequences.ensembl_gene_id|hsapiens_gene_ensembl.default.sequences.ensembl_transcript_id|hsapiens_gene_ensembl.default.sequences.3utr&FILTERS=&VISIBLEPANEL=resultspanel

# <headingcell level=5>

# Collect UTR sequences

# <codecell>

# format query so that it can be used via biomart REST service
# http://www.biomart.org/martservice.html

try:
    ## request data via HTTP REST service
    queryurl = baseurl + '?' + urlencode({'query': query})
    filename = wget.download(queryurl, out=OUTPUT_DIR + 'hsap.3utr.fasta', bar="")
except:
    # SSL problems means that this may not work on windows:
    # requires this workaround - 
    # http://stackoverflow.com/questions/13167907/https-request-results-in-reset-connection-in-windows-with-python-3
    # related to this bug: http://bugs.python.org/issue16361
    
    print("wget failed, trying urllib method")
    ssl_context = urllib.request.HTTPSHandler(
                                              context=ssl.SSLContext(ssl.PROTOCOL_TLSv1))
    opener = urllib.request.build_opener(ssl_context)
    urllib.request.install_opener(opener)
    queryurl = Request(baseurl, urlencode({'query': query}).encode('ascii'))
    
    # do request / download file
    # queryfile = wget.download(urlopen(queryurl, timeout=1000), out=OUTPUT_DIR, bar="")
    filename = OUTPUT_DIR + 'hsap.3utr.fasta'
    with open(filename, mode='wt', ) as f:
        g = urlopen(queryurl, timeout=1000)
        f.write(g.read())

# <codecell>

try:
    print(filename)
except:
    filename = OUTPUT_DIR + 'hsap.3utr.fasta'
    print(filename)

# <headingcell level=5>

# Parser for UTR FASTA file

# <codecell>

# adapted from: http://stackoverflow.com/questions/7654971/parsing-a-fasta-file-using-a-generator-python
def read_fasta(fp, limit=None):
    ### return tuple:
    ### id line : sequence
    
    # initialise holders
    name, seq = None, []
    count = 0
    # open file and go through one line at a time
    for line in fp:
        # finish if limit reached (number of records)
        if limit and count >= limit:
            break
        line = line.rstrip()
        if line.startswith('Seq'): # "Sequence not available"
           # print(line)
            continue
        # if we find a new ID line
        if line.startswith(">"):
            count = count + 1 # record count of new id lines
            # yeild results so far (previous record if it exists)
            if name and len(''.join(seq)) > 0:
                yield (name, ''.join(seq))
            # reset holders with new id
            name, seq = line, []
        else:
            # otherwise add sequence line to list
            seq.append(line)
    # return the last record
    if name: 
        yield (name, ''.join(seq))

# <headingcell level=5>

# Split UTR seqences to create sequence : UTR node pairs

# <codecell>

# read and parse entries
utr = {}
# collect id : sequence pairs
with open(filename) as f:
    for name, seq in read_fasta(f, limit =5e3): # 1e5 = ~1/2 UTRs
        utr[name] = seq 
print('' + str(len(utr.keys())) + ' UTR sequences read in for parsing')

# parse full sequences to N-mers
for key in utr.keys():      
    # split to individual characters
    # replace single characters with 'sliding window' N-mer
    tmp = {
           'l':len(utr[key]), # length of utr
           't':'utr', # type
           's':[], # list of seeds
           'gene': key.strip('>').split('|')[0], # gene id
           'transcript': key.strip('>').split('|')[1] # transcript id
           }
    for i in range(0, len(utr[key])-(NMER_LEN)):
        tmp['s'].append( utr[key][i:i+(NMER_LEN)] )
    utr[key] = tmp

if True:
    # check a few keys
    for i, k in enumerate(utr.keys()):
        if i >= 2: break
        print(k)
        for k2, v2 in utr[k].items():
            if k2 !='s' :
                print(k2 + ':' + str(v2))
        print(len(utr[k]['s']))
        print(utr[k]['s'][:3])

# <headingcell level=2>

# Create Graph object

# <headingcell level=5>

# Load UTR : sequence pairs to networkx object

# <rawcell>

# Aim is to produce a miltipartite graph consisting of:
# 
# Gene --> transcript  - (3'UTR) -> seed
# 
# i.e. :
# Gene1 -> transcript1 -> seed1
# Gene1 -> transcript1 -> seed2
# Gene1 -> transcript2 -> seed1
# Gene1 -> transcript2 -> seed2
# Gene1 -> transcript2 -> seed3
# Gene2 -> transcript1 -> seed1
# ...
# 
# _Gene_
#  - node id = ensembl gene id
#  - type = gene
#  - size = 10
# _Transcript_ 
#  - node id = ensembl transcript id
#  - type = transcript
#  - size = 5
# - other attr: length of UTR sequence (for filtering?)
# _seed_ 
#  - node id = seed sequence
#  - type = seed
#  - size = 2
# 
# NB: does each ensembl transcript id have only 1 utr? i.e. if the utr changes do you get a new transript id?

# <codecell>

# initiate graph object (undirected, no multi-edges)
G = nx.Graph()

# <codecell>

# collect node ids in dict for intitialisation and lookups
nodes = {
         # get unique set of genes and transcripts
         'genes': set([ utr[k]['gene'] for k in utr.keys()]), 
         'transcripts' : set([ utr[k]['transcript'] for k in utr.keys()]), 
         # flatten 2d list and convert to set for unique seed sequences
         'seeds' : set([element for sub in [utr[k]['s'] for k in utr.keys()] for element in sub]) 
         }

# check numbers
for k in nodes.keys():
    print(str(len(nodes[k])) + ' node ids collected for ' + k)

# <codecell>

# load gene nodes as
# node, attribute dict tuples list
G.add_nodes_from([(k, {'type' : 'gene', 'size' : 200, '_id' : k, 'label' : k, 'col' : 'w'}) for k in nodes['genes']])

# <codecell>

# load UTR (transcript) nodes
G.add_nodes_from([(k, {'type' : 'transcript', 'size' : 100, '_id' : k, 'label' : k, 'col' : 'b'}) for k in nodes['transcripts']]) 

# <codecell>

# load seed nodes
G.add_nodes_from([(k, {'type' : 'seed', 'size' : 25, '_id' : k, 'label' : '', 'col' : 'c'}) for k in nodes['seeds']]) 

# <codecell>

## clean up
# delete nodes lists
del nodes

# <codecell>

# add edges gene --> transcript
G.add_edges_from([(utr[k]['gene'], utr[k]['transcript']) for k in utr.keys()])

# <codecell>

# add edges transcript --> seed
for i, v in enumerate(utr.values()):
    G.add_edges_from([(v['transcript'], seed) for seed in v['s']])
    v = {}
    if i % 5e3 == 0:
        print(str(i) + ' UTRs with all seed interactions added to graph')

# <codecell>

## clean up
# remove utr
del utr

# <headingcell level=2>

# Analysis & visualisation

# <codecell>

# add degree as node attribute
nx.set_node_attributes(G,'degree',{key:value for key, value in G.degree().items()})

# <codecell>

## plot degree distribution by node type: can we prune 'leaf' seed nodes to reduce graph size?
# collect data
tmp = {}
for n in G.node.values():
    if n['type'] not in tmp:
        tmp[n['type']] = []
    tmp[n['type']].append( n['degree'] )
    
## plot as box/whisker
# initialise plot & axes
fig, ax = plt.subplots(1, len(tmp.keys()), figsize=(18,18))
for i, k in enumerate(tmp.keys()):
    # Create the boxplot
    ax[i].boxplot(tmp[k], 1)
    # x-axis labels
    ax[i].set_xticklabels(k)

plt.show()

# <codecell>

def generic_graph_plot(G):
    ## simple plot of all nodes
    # if graph is small enough
    if len(G.node.values()) < 4000:
        # open figure
        plt.figure(figsize=(18,18))
        
        # get parameters for plotting
        cols = list(i['col'] for i in G.node.values())
        size = list(i['size'] for i in G.node.values())
        labels = {i['_id']:i['label'] for i in G.node.values()}
        
        nx.draw_networkx(G,
                pos=nx.spring_layout(G, weight=1, scale=10.0, iterations=50),
                node_color=cols, node_size=size,
                with_labels=False, alpha=0.5, labels=labels
                )
    else:
        print('Graph too large to print it all: consider downsampling/subsetting/derivations')

# <codecell>

def weighted_graph_plot(G, threshold_property = 'weight', threshold = 0, scale=50, debug=False):
    
    ### plot a weighted layout of graph <G>
    # parameter G: input networkx graph object
    # parameter threshold_property: edge property on which thresholding is to be done (will also be used as edge weight for plotting)
    # parameter threshold: int, minimum edge weight value to be considered
    # parameter scale: int, for plotting scale all edge thicknesses (weight) to this value
    # parameter debug: default false, print debugging info?
    
    import numpy as np
    
    ## simple plot of all nodes
    
    # open figure
    plt.figure(figsize=(18,18))
    
    # create subgraph only containing edges with weight above threshold
    # and nodes with a degree of at least 1 (i.e. no disconnected nodes after removing edges below threshold)
    H = G.copy()
    
    if debug:
        # print out some debugging info
        print('total starting nodes: ' + str(len(H.nodes())))
        print('total starting edges: ' + str(len(H.edges())))
        
        print('edges to remove: ' + str(len([ (i[:2]) for i in H.edges(data=True) if i[2][threshold_property] < threshold ])))
    
    H.remove_edges_from([ (i[:2]) for i in H.edges(data=True) if i[2][threshold_property] < threshold ])
    H.remove_nodes_from([ k for k,v in H.degree().items() if v ==0 ])
    
    if debug:
        # print out some debugging info
        print('total nodes after thresholding: ' + str(len(H.nodes())))
        print('total nodes after thresholding with 0 degree: ' + str(len([ k for k,v in H.degree().items() if v ==0 ])))
        print('total edges after thresholding: ' + str(len(H.edges())))
        
        print('Nodes:')
        print(H.nodes(data=True))
        print('Edges:')
        print(H.edges(data=True))
    
    # if graph is small enough
    # parameter threshold: minimum edge weight for plotting (assumed to be in edge property called 'weight') 
    if len(H.nodes()) < 2000:
    
        # get parameters for plotting
        cols = list(i['col'] for i in H.node.values())
        size = list(i['size'] for i in H.node.values())
        labels = {i['_id']:i['label'] for i in H.node.values()}
        width = list(i[2][threshold_property] for i in H.edges(data=True))
        
        # scale width for plotting
        width = [ (i*scale)/np.amax(width) for i in width ]
        
        nx.draw_networkx(H,
                pos=nx.spring_layout(H, weight='weight', scale=10.0, iterations=50), # Fruchterman-Reingold force-directed
                node_color=cols, node_size=size, width=width,
                # only plot nodes present in edgelist
                # nodelist=list(set( flatten([ (i[:2]) for i in G.edges(data=True) if i[2]['weight'] >= threshold ]) )),
                # only plot edges with weight of at least threshold
                # edgelist=[ (i[:2]) for i in G.edges(data=True) if i[2]['weight'] >= threshold ],
                with_labels=False, alpha=0.5, labels=labels
                )
        
    else:
        print('Graph too large to print it at ' + threshold_property + ' > ' + threshold + ' : consider downsampling/subsetting/derivations')

# <codecell>

generic_graph_plot(G)

# <codecell>

# common neighbours function for cocitation
def common_neighbours(graph, n1, n2):
    ### calculate the common neighbours for nodes <n1> and <n2> in graph <graph>
    ### return list of common neighbours, or empty list if none
    
    # NB: assumes n1 and n2 are in <graph>, DOES NOT CHECK
    cn = set()
    if n1!=n2: # no self loops
        cn = set(graph.neighbors(n1)) & set(graph.neighbors(n2))
    
    return(list(cn))

# <codecell>

# common neighbours function for cocitation
def all_neighbours(graph, nbunch):
    ### calculate all the neighbours for nodes in <nbunch> in graph <graph>
    ### return list of neighbours, excluding <nbunch>, or empty list if none
    
    # NB: assumes nbunch are in <graph>, DOES NOT CHECK
    
    an = list()
    for n in nbunch:
        an.append([ graph.neighbors(n) ])
    an = set(flatten(an)) - set(nbunch) # remove nbunch from neighbour list
    
    return(list(an))

# <codecell>

# function for cocitation between two nodes in graph, using helper functions <all_neighbours()> and <common_neighbours()>
def cocite(x):
    (graph, n1, n2) = x
    weight = len(common_neighbours(graph, n1, n2))
    
    if weight > 0:
        # include jaccard as a scaled measure : http://en.wikipedia.org/wiki/Jaccard_index
        jaccard_index = (weight / len(all_neighbours(graph, [n1, n2])))
        jaccard_dist = 1- jaccard_index
        
        # return tuple for adding cocitation edge
        return(tuple([n1, n2, {'weight' : weight, 'j_dist' : jaccard_dist, 'j_idx' : jaccard_index}]))
    
    else:
        return(None)

# <codecell>

# cocitation wrapper function
def cocitation(graph, nbunch=None, debug=False, matrix=False, parallel=None):
    ### for nodes in <nbunch> in graph <graph>
    ### calculate cocitation graph, edges weighted by number of common nodes
    ### return new networkx graph
    
    # parameter graph: must be a networkx graph instance
    # parameter nbunch: must be node ids from graph <graph>
    #                   node ids not in graph <graph>, then they will be ignored
    #                   default = None (calculate for all nodes in graph)
    # parameter debug: default false, print debugging info?
    # parameter matrix: return an adjacency matrix? if false (default), returns networkx.Graph object
    
    # TODO: param checking and error checking
    
    import networkx as nx
    import numpy as np
    
    if debug:
        print('Checking input nbunch list')
    # if no nbunch supplied, do for all nodes
    if nbunch is None:
        print('Calculating cocitation for all nodes')
        nbunch = [k for k in graph.node.keys()]
    
    # drop any nodes not in graph
    nbunch = list(set(nbunch) & set(graph.node.keys())) # set intersection
    
    if debug:
        # check input list
        print('nbunch type: ' + str(type(nbunch)))
        print('nbunch length: ' + str(len(nbunch)))
        print(nbunch[:5])
        if len(set(nbunch)) != len(nbunch):
            print('nbunch is not unique!')
    
    # new graph object to be returned, copying over node attributes from old graph
    if debug:
        print('Creating cocitation graph')
    H = nx.Graph()
    H.add_nodes_from([n for n in graph.nodes(data=True) if n[0] in nbunch])
    
    if debug:
        # check graph
        print('number of nodes in input graph: ' + str(len(graph.nodes())))
        print('number of nodes in input nbunch: ' + str(len(nbunch)))
        print('number of nodes in output cocitation graph: ' + str(len(H.nodes())))
    
    if debug:
        print('Adding cocitation edges to graph')
    if not parallel:
        # iterate through all nodes in nbunch
        # adding an edge if they are cocited, with edge weight as a 'weight' attribute
        # do cocitation
        e = map(cocite, 
                # only cover 'lower half' of adjacency matrix
                [(graph, nbunch[i], nbunch[j]) for i in range(len(nbunch)) for j in range(len(nbunch)) if i < j ]
                )
        H.add_edges_from([i for i in e if i]) # if no cocitation, cocite() return None
                
    else:
        # iterate through node pairs as above, but in parallel
        from multiprocessing import Pool
        
        # init processes
        p = Pool(8)
        
        e = p.map(cocite, [(graph, nbunch[i], nbunch[j]) for i in range(len(nbunch)) for j in range(len(nbunch)) if i < j ])
        
        H.add_edges_from([i for i in e if i]) # if no cocitation, cocite() return None
    
    if matrix:
        return(to_numpy_matrix(H, nodelist=nbunch, dtype=int, weight='weight', nonedge=0))
    else:
        return(H)

# <codecell>

# test
generic_graph_plot(
                   cocitation(G, nbunch=['ENST00000225504', 'ENST00000220509', 'ENST00000238508', 'ENST00000296849', 'ENST00000296847'], 
                              debug=True)
                   )

# <codecell>

# try it out
cocit = cocitation(G, nbunch=[k for k,v in G.node.items() if v['type'] == 'transcript'], debug=True) # very slow & wasteful
# cocit = cocitation(G, nbunch=[k for k,v in G.node.items() if v['type'] == 'transcript'], debug=True, parallel=True)

# <codecell>

# save cocitation for later use
nx.write_gml(cocit, OUTPUT_DIR + 'hsap.3utr.cocit.gml.gz')
# cocit = nx.read_gml(OUTPUT_DIR + 'hsap.3utr.cocit.gml', relabel=True)

# <codecell>

## plot as scatter: jaccard vs weight
# initialise plot & axes
fig, ax = plt.subplots(1, 1, figsize=(8,8))
# create plot
ax.scatter([i[2]['j_dist'] for i in cocit.edges(data=True)], [i[2]['weight'] for i in cocit.edges(data=True)],
           alpha=.5)
ax.set_ylabel('edge weight')
ax.set_xlabel('Jaccard distance')
plt.show()

# <codecell>

# plot as graph, varying weight threshold
weighted_graph_plot(cocit, threshold=0)
weighted_graph_plot(cocit, threshold=10)
weighted_graph_plot(cocit, threshold=100, debug=True)

