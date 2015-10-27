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
    for name, seq in read_fasta(f, limit=2e4): # ~80k UTRs total
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
    print('Dictionary structure after parsing:')
    for i, k in enumerate(utr.keys()):
        if i >= 2: break
        print(k)
        for k2, v2 in utr[k].items():
            if k2 !='s' :
                print(k2 + ':' + str(v2))
        print(len(utr[k]['s']))
        print(utr[k]['s'][:3])

# <codecell>

# invert dictionary to create {seed : [list of transcripts], ...} structure
# i.e. ready for creation of cocitation graph, without intermediate bipartite graph a la biomart_test

seed = dict(
            # flatten 2d list and convert to set for unique seed sequences
             [(element, []) for sub in [utr[k]['s'] for k in utr.keys()] for element in sub]
             )

for v in utr.values():
    for s in v['s']:
        seed[s].append(v['transcript'])

# in place sort of list - important for adding edges later
for v in seed.values():
    v.sort()

# <headingcell level=2>

# Create Graph object

# <headingcell level=5>

# Construct UTR cocitation network

# <rawcell>

# Aim is to produce a miltipartite graph consisting of:
# 
# Gene --> transcript  <- (cocitation (i.e. common seeds) -> transcript
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
# becomes :
# 
# G1 -> t1 -> t2 -> G1
# G2 -> t1 -> t1 -> G1
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
# 
# NB: does each ensembl transcript id have only 1 utr? i.e. if the utr changes do you get a new transript id?

# <codecell>

## initiate graph object
# undirected, no multi-edges
# latter point is important in adding edges to graph object
# in nx.Graph() node pair can only be connected once
# in nx.MultiGraph() they can have multiple edges
G = nx.Graph()

# <codecell>

# collect node ids in dict for intitialisation and lookups
nodes = {
         # get unique set of genes and transcripts
         'genes': set([ utr[k]['gene'] for k in utr.keys()]), 
         'transcripts' : set([ utr[k]['transcript'] for k in utr.keys()])
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

## clean up
# delete nodes lists
del nodes

# <codecell>

# add edges gene --> transcript
G.add_edges_from([(utr[k]['gene'], utr[k]['transcript'], {'type':'gene_transcript'}) for k in utr.keys()])

# <headingcell level=5>

# Calculate cocitation for all UTRs

# <codecell>

# add edges transcript --> seed
print('UTR -> UTR cocitation edges being added to graph, seeds being processed:')
for i, s in enumerate(seed.values()):
    G.add_edges_from([(s[i], s[j], {'weight': 0, 'type': 'cocitation'}) for i in range(len(s)) for j in range(len(s)) if i < j])
    if i % abs((4^NMER_LEN)/10) == 0: # ~10 updates
        print(str(i) + ' seeds processed')
print(str(len(G.edges())) + ' edges added to graph')
print(G.edges(data=True)[:3])

# <codecell>

# reiterate through seed -> transcripts, updating edge weight
print('UTR -> UTR cocitation edges being updated, seeds being processed:')
for i, s in enumerate(seed.values()):
    for (x,y) in [(s[i], s[j]) for i in range(len(s)) for j in range(len(s)) if i < j]:
        G[x][y]['weight'] = G[x][y]['weight'] + 1
    if i % abs((4^NMER_LEN)/10) == 0:
        print(str(i) + ' seeds processed')
print(str(len(G.edges())) + ' edges updated')
print(G.edges(data=True)[:3])

# <codecell>

## one final cycle to update cocitation edges with jaccard stats
# translation table for UTR : transcripts for lookups
tmp = dict([(v['transcript'], k) for k,v in utr.items()]) 

# iterate through and update edge properties
print('UTR -> UTR cocitation edges being updated:')
for n, (i,j) in enumerate(G.edges()):
    if G[i][j]['type'] == 'cocitation':
        # only deal with cocitation edges
        
        G[i][j]['j_idx'] = G[i][j]['weight'] / len(set(utr[tmp[i]]['s']).union(set(utr[tmp[j]]['s'])))
        G[i][j]['j_dist'] = 1 - G[i][j]['j_idx']
        
    else:
        continue
    
    if n % 1e5 == 0:
        print(str(n) + ' edges updated')

# clean up
del tmp
# check results
print(G.edges(data=True)[:5])

# <codecell>

# check results
tmp = dict([(v['transcript'], k) for k,v in utr.items()])
for (i,j,d) in G.edges(data=True):
    
    if d['type'] == 'cocitation' and d['j_dist'] < 0.2 and d['weight'] < 100:
        print(d)
        print(utr[tmp[i]])
        print(utr[tmp[j]])
del tmp

# <headingcell level=2>

# Analysis & visualisation

# <codecell>

# add degree as node attribute
nx.set_node_attributes(G,'degree',{key:value for key, value in G.degree().items()})
nx.set_node_attributes(G,'weighted_degree',{key:value for key, value in G.degree(weight='weight').items()})

# <codecell>

def generic_graph_plot(G):
    ## simple plot of all nodes
    # if graph is small enough
    if len(G.node.values()) < 4000:
        # open figure
        p = plt.figure(figsize=(18,18))
        
        # get parameters for plotting
        cols = list(i['col'] for i in G.node.values())
        size = list(i['size'] for i in G.node.values())
        labels = {i['_id']:i['label'] for i in G.node.values()}
        
        nx.draw_networkx(G,
                pos=nx.spring_layout(G, weight=1, scale=10.0, iterations=50),
                node_color=cols, node_size=size,
                with_labels=False, alpha=0.5, labels=labels
                )
        return(p)
    else:
        print('Graph too large to print it all: consider downsampling/subsetting/derivations')
        return(None)

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
    # if graph is small enough
    # parameter threshold: minimum edge weight for plotting (assumed to be in edge property called 'weight') 
    if len(G.node.values()) < 4000:
        # open figure
        p = plt.figure(figsize=(18,18))
        
        # create subgraph only containing edges with weight above threshold
        # and nodes with a degree of at least 1 (i.e. no disconnected nodes after removing edges below threshold)
        H = G.copy()
        
        if debug:
            # print out some debugging info
            print('total starting nodes: ' + str(len(H.nodes())))
            print('total starting edges: ' + str(len(H.edges())))
            
            print('edges to remove: ' + str(len([ (i[:2]) for i in H.edges(data=True) if i[2].get(threshold_property, -1) < threshold ])))
        
        H.remove_edges_from([ (i[:2]) for i in H.edges(data=True) if i[2].get(threshold_property, -1) < threshold ])
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
        return(p)
    else:
        print('Graph too large to print it all: consider downsampling/subsetting/derivations')
        return(None)

# <codecell>

generic_graph_plot(G)
weighted_graph_plot(G, threshold_property='weight', threshold=100)
weighted_graph_plot(G, threshold_property='j_idx', threshold=0.5)

# <codecell>

## plot as scatter: jaccard vs weight
# initialise plot & axes
fig, ax = plt.subplots(1, 1, figsize=(8,8))
# create plot
ax.scatter([d['j_idx'] for (i,j,d) in G.edges(data=True) if d['type'] == 'cocitation'], 
           [d['weight'] for (i,j,d) in G.edges(data=True) if d['type'] == 'cocitation'])
ax.set_xlim([0,1])

plt.show()
plt.savefig(OUTPUT_DIR + 'JI_vs_weight.png')

# <codecell>

for thresh in [i/100 for i in range(10,35, 5)]:
    p = weighted_graph_plot(G, threshold_property='j_idx', threshold=thresh)
    p.suptitle('Jaccard index threshold: ' + str(thresh), fontsize=25)
    p.savefig(OUTPUT_DIR + 'cocit.JI_gt_' + str(thresh) + '_.png')
    

# <codecell>

# save cocitation for later use
nx.write_gml(G, OUTPUT_DIR + 'hsap.3utr.cocit.gml.gz') # automajickally compressed
# G = nx.read_gml(OUTPUT_DIR + 'hsap.3utr.cocit.gml', relabel=True)

