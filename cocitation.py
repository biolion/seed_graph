# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # UTR 'cocitation' networks
# Using the number of shared N-mer sequences
# construct a graph linking UTRs (and genes)
# 
#  - how similar (sequence-wise) are UTRs?
#  - Do we see any strongly similar UTRs?
#  - Do they share biological properties?
#  - Are they common 'off-target partners'?

# <markdowncell>

# #### About
# _authors_
# 
#  - clayie1
#  - nigscfl1
# 
# _version_
# 
#  - August 2014
#  - commit d1fe81e89da8c3712b804678a176c18c36582c2a

# <markdowncell>

# ### Set up environment

# <codecell>

# general
import os
import sys
import urllib.request
from urllib.request import Request, urlopen
from urllib.parse import urlencode
import wget
# import ssl # required for request to biomart on windows

# <codecell>

# analysis
import pandas as pd
import networkx as nx
import numpy as np
import scipy as sp
import scipy.spatial as spspat
import scipy.cluster as spclu
import scipy.sparse as sps

# <codecell>

# plotting
import matplotlib as mpl
import matplotlib.pyplot as plt

# <codecell>

## OPTION FLAGS

# how long should Nmers be when parsing UTR sequences
NMER_LEN = 8

# Run on test dataset?
TEST = False

# min/max UTR length (None = no limit)
UTR_RANGE = (NMER_LEN + 9, 10000) # min: must be able to produce 10 seeds

# min/max seed representation (None = no limit)
SEED_RANGE = (None, 5000)

# use cached downloaded UTR file?
USE_UTR_CACHE = True
if(USE_UTR_CACHE):
    print('Warning: Fresh data will not be downloaded, using cached data')
    
# use cached cocitation files?
USE_COCIT_CACHE = True
if(USE_COCIT_CACHE):
    print('Warning: Cocitation will not be recalculated, using cached data')

# print debugging / checks?
VERBOSE = True
LOGGING = False
if (LOGGING):
    sys.stdout = open(OUTPUT_DIR + 'logfile.txt', "w")

# <codecell>

# global vars
DIR = os.getcwd()
OUTPUT_DIR = DIR + '/' + 'output/'
TMP_DIR = '/da/dmp/cb/clayie1/tmp/seed_graph/'
if not os.path.exists(OUTPUT_DIR):
    print('Making ' + OUTPUT_DIR)
    os.makedirs(OUTPUT_DIR)
print('Output in ' + OUTPUT_DIR)
# PROXY = os.environ.get('http_proxy', 'http://eu-chbs-PROXY.eu.novartis.net:2011/')
PROXY = os.environ.get('http_proxy', 'http://nibr-proxy.global.nibr.novartis.net:2011/')

# <codecell>

# parallel support
from IPython.parallel import Client, error
cluster = Client() # default profile
dview = cluster.direct_view() # direct access/control (including push/pull)
lbview = cluster.load_balanced_view() # load balanced view for running jobs
cluster.ids
# nb: MPI clusters must be initiated - see: http://nbviewer.ipython.org/github/ipython/ipython/blob/master/examples/Parallel%20Computing/Using%20MPI%20with%20IPython%20Parallel.ipynb

# test
dview.execute('import os')
%px print("testing... pid: " + str(os.getpid()))

# <markdowncell>

# ## Collect UTR sequences from Biomart

# <markdowncell>

# ### Set up query

# <codecell>

if not (USE_UTR_CACHE):
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

# <markdowncell>

# ### Collect UTR sequences

# <codecell>

# destination of downloaded data / cached data
filename = OUTPUT_DIR + 'hsap.3utr.fasta'
if not (USE_UTR_CACHE):
    # format query so that it can be used via biomart REST service
    # http://www.biomart.org/martservice.html

    try:
        ## request data via HTTP REST service
        queryurl = baseurl + '?' + urlencode({'query': query})
        filename = wget.download(queryurl, out=filename, bar="")
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
        with open(filename, mode='wt', ) as f:
            g = urlopen(queryurl, timeout=1000)
            f.write(g.read())

# <markdowncell>

# ### Parser for UTR FASTA file

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

# <markdowncell>

# ### Collect and parse UTR sequences

# <codecell>

# read and parse entries
utr = {}

# collect id : sequence pairs
with open(filename) as f:
    for name, seq in read_fasta(f, limit=[None, 1e3][TEST + 0]): # ~80k UTRs total, if TEST read them all
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

# <codecell>

# check a few keys
if VERBOSE:
    print('Dictionary structure after parsing:')
    for i, k in enumerate(utr.keys()):
        if i >= 2: break
        print(k)
        for k2, v2 in utr[k].items():
            if k2 !='s' :
                print(k2 + ':' + str(v2))
        print(len(utr[k]['s']))
        print(utr[k]['s'][:3])

# <markdowncell>

#     what about redundancy of seeds? i.e. if UTR A contains the seed 'AATGAATG' 12 times and
#     UTR B contains the same seed 3 times, is that a 'seed match' of 1, 3, 12?

# <markdowncell>

# ### Examine and filter UTRs based on size (if given in FLAGS)

# <codecell>

# check size distribution of UTRs
# initialise plot
p = plt.figure(figsize=(12,8))
# create plot
plt.hist([utr[k]['l'] for k in utr.keys()], 100, facecolor='green', alpha=0.5, normed=True)
counts, bin_edges = np.histogram([utr[k]['l'] for k in utr.keys()], 100, normed=True)
cdf = np.cumsum(counts)
cdf = (cdf / np.max(cdf)) * np.max(counts) # scale to biggest bin
plt.plot(bin_edges[1:], cdf)
plt.xlabel('UTR length')
plt.ylabel('Probability')
plt.suptitle('UTR length distribution')
if VERBOSE:
    plt.show()
p.savefig(OUTPUT_DIR + 'UTR_length_distribution.png')

# <codecell>

# filter UTRs
if (UTR_RANGE[0] is not None or UTR_RANGE[1] is not None):
    # i.e. only do filtering if one or both range limts are set
    if (UTR_RANGE[0] is None):
        UTR_RANGE = (0, UTR_RANGE[1]) # i.e. no lower limit
    if (UTR_RANGE[1] is None):
        UTR_RANGE = (UTR_RANGE[0], np.max([v['l'] for v in utr.values()])) # set to max, i.e. no limit
    
    # do filter
    utr = dict((k,v)for (k,v) in utr.items() if UTR_RANGE[0] <= v['l'] <= UTR_RANGE[1])
if VERBOSE:
    print(str(len(utr.keys())) + ' UTRs remaining after filtering')

# <markdowncell>

# ### Split UTR seqences to create sequence : UTR node pairs

# <codecell>

# invert dictionary to create {seed : [list of transcripts], ...} structure
# i.e. ready for creation of cocitation graph, without intermediate bipartite graph a la biomart_test
print('Creating seeds from parsed UTRs')
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
if VERBOSE:
    print(str(len(seed.keys())) + ' seeds extracted')

# <markdowncell>

# ### examine seed representation and filter

# <codecell>

# distribution of transcripts per seed
# which seeds are present in the most transcripts?
# initialise plot
p = plt.figure(figsize=(12,8))
# create plot
plt.hist([len(v) for v in seed.values()], 100, facecolor='green', alpha=0.5, normed=True, log=True)
counts, bin_edges = np.histogram([len(v) for v in seed.values()], 100, normed=True)
cdf = np.cumsum(counts)
cdf = (cdf / np.max(cdf)) * np.max(counts) # scale to biggest bin
plt.plot(bin_edges[1:], cdf)
plt.xlabel('Transcripts per seed')
plt.ylabel('Probability')
plt.suptitle('Seed representation distribution')
if VERBOSE:
    plt.show()
p.savefig(OUTPUT_DIR + 'seed_representation_distribution.png')

# <codecell>

# table of most common seeds
df = pd.DataFrame.from_items([('seed', [k for k in seed.keys()]), ('representation', [len(v) for v in seed.values()])])
df.sort('representation', ascending=False, inplace=True)
if VERBOSE:
    print(df.head(n=20))
df.to_csv(OUTPUT_DIR + str(NMER_LEN) + 'mer_seed_representation.csv', index=False)

# <codecell>

# filter seeds
if (SEED_RANGE[0] is not None or SEED_RANGE[1] is not None):
    # i.e. only do filtering if one or both range limts are set
    if (SEED_RANGE[0] is None):
        SEED_RANGE = (0, SEED_RANGE[1]) # i.e. no lower limit
    if (SEED_RANGE[1] is None):
        SEED_RANGE = (SEED_RANGE[0], np.max([len(v) for v in seed.values()])) # set to max, i.e. no limit
    
    # do filter
    seed = dict((k,v)for (k,v) in seed.items() if SEED_RANGE[0] <= len(v) <= SEED_RANGE[1])
if VERBOSE:
    print(str(len(seed.keys())) + ' seeds remaining after filtering')

# <markdowncell>

# ## Calculate cocitation for all UTRs

# <markdowncell>

# ### prepare for parallel evaluation of transcript-transcript cocitation

# <codecell>

## clear/create temp directory for storing intermediate results
print('UTR -> UTR cocitation edges being calculated:')

# make/clear tmp directory
if not os.path.exists(TMP_DIR):
    if USE_COCIT_CACHE:
        print('Warning: Cocitation caching enabled, but cache not present. Overriding flag and calculating fresh.')
        USE_COCIT_CACHE = False
    os.makedirs(TMP_DIR)
else:
    # remove existing/old  temp files if we are not caching
    if not USE_COCIT_CACHE:
        for f in [f for f in os.listdir(TMP_DIR) if f.endswith('.tmp')]:
            os.remove(os.path.join(TMP_DIR, f))

if not USE_COCIT_CACHE:
    ## prepare data needed by each process
    # see: http://stackoverflow.com/questions/11371009/parallel-mapping-functions-in-ipython-w-multiple-parameters
    print('Preparing parallel environment')
    # required imports
    # dview.execute('import os')
    # dview.execute('import panadas as pd')
    # transcripts to be processed
    todo = [k for k in utr.keys()] # all processes get a todolist
    # convert to tuple (immutable!)
    todo = tuple(todo)
    
    # push data to processes using direct view
    dview.push(dict(utr=utr, todo=todo, TMP_DIR=TMP_DIR, NMER_LEN=NMER_LEN))

# <markdowncell>

# ### parallel function
# http://ipython.org/ipython-doc/stable/parallel/parallel_task.html

# <codecell>

@lbview.parallel(block=True)
def cocit(i):
    
    import os
    import pandas as pd
    
    print(' > starting cocitation for %s, (pid: %i)', todo[i], os.getpid())
    
    df = list()
    
    # start cocitation analysis
    for j in range(i+1, len(todo), 1):
        
        # weight (# common seeds), jaccard index and distance
        weight = len( set(utr[todo[i]]['s']).intersection(set(utr[todo[j]]['s'])) )
        if weight >= NMER_LEN * 2: ### need to edit this...
            j_idx = weight / len(set(utr[todo[i]]['s']).union(set(utr[todo[j]]['s'])))
            j_dist = 1 - j_idx
        
            # add to df
            # to, from, weight, ji, jd
        
            df.append([utr[todo[i]]['transcript'], utr[todo[j]]['transcript'], weight, j_idx, j_dist])
        
    # combine and print results
    df = pd.DataFrame(df, columns=['to', 'from', 'weight', 'j_idx', 'j_dist'])
    df.to_csv(os.path.join(TMP_DIR, utr[todo[i]]['transcript'] + '.tmp'), index=False)
    
    return(i)

# <codecell>

# run it
if not USE_COCIT_CACHE:
    print(str(np.max(cocit.map(range(len(todo) -1)))) + ' transcripts processed')
else:
    print('Using cached cocitation data')

# <markdowncell>

# # Analysis & visualisation

# <markdowncell>

# ### read back in .tmp files and collect results

# <codecell>

results = None
for f in [os.path.join(TMP_DIR, f) for f in os.listdir(TMP_DIR) if f.endswith('.tmp')]:
    tmp = pd.DataFrame.from_csv(f, index_col=[0,1])
    tmp = tmp[tmp['weight'] > 0] # remove non edges
    if results is None:
        # first file, so start the results data frame
        results = tmp.copy()
    else:
        # combine results
        results = results.append(tmp, ignore_index=False)

# check results
print(results.shape)

results.head()

# <markdowncell>

# ### plotting jaccard and weight statistics

# <codecell>

# check size distribution of weight
# initialise plot
p = plt.figure(figsize=(12,8))
# create plot
plt.hist(results['weight'], 100, facecolor='green', alpha=0.5, normed=True)
counts, bin_edges = np.histogram(results['weight'], 100, normed=True)
cdf = np.cumsum(counts)
cdf = (cdf / np.max(cdf)) * np.max(counts) # scale to biggest bin
plt.plot(bin_edges[1:], cdf)
plt.xlabel('Weight')
plt.ylabel('Probability')
plt.suptitle('Weight distribution')
if VERBOSE:
    plt.show()
p.savefig(OUTPUT_DIR + str(NMER_LEN) + 'mer_weight_distribution.png')

# <codecell>

# check size distribution of weight
# initialise plot
p = plt.figure(figsize=(12,8))
# create plot
plt.hist(results['j_idx'], 100, facecolor='green', alpha=0.5, normed=True)
counts, bin_edges = np.histogram(results['j_idx'], 100, normed=True)
cdf = np.cumsum(counts)
cdf = (cdf / np.max(cdf)) * np.max(counts) # scale to biggest bin
plt.plot(bin_edges[1:], cdf)
plt.xlabel('Jaccard Index')
plt.ylabel('Probability')
plt.suptitle('Jaccard Index distribution')
if VERBOSE:
    plt.show()
p.savefig(OUTPUT_DIR + str(NMER_LEN) + 'mer_JI_distribution.png')

# <codecell>

## as scatter
# initialise plot & axes
fig, ax = plt.subplots(1, 1, figsize=(12,12))
# create plot
ax.scatter(results['j_idx'], results['weight'])
ax.set_xlim([0,1])
ax.set_ylim(bottom=0)
ax.set_xlabel('Jaccard Index')
ax.set_ylabel('Weight')
if VERBOSE:
    plt.show()
fig.savefig(OUTPUT_DIR + str(NMER_LEN) + 'mer_JI_vs_weight.png')

# <markdowncell>

# ## Construct UTR cocitation network

# <markdowncell>

#     Aim is to produce a *miltipartite graph* consisting of:
#     _Gene --> transcript  <- (cocitation (i.e. common seeds) -> transcript_
#         i.e. :
#             Gene1 -> transcript1 -> seed1
#             Gene1 -> transcript1 -> seed2
#             Gene1 -> transcript2 -> seed1
#             Gene1 -> transcript2 -> seed2
#             Gene1 -> transcript2 -> seed3
#             Gene2 -> transcript1 -> seed1
#             ...
#         becomes :
#             G1 -> t1 -> t2 -> G1
#             G2 -> t1 -> t1 -> G1
# 
#         *Nodes*:
#             _Gene_
#             - node id = ensembl gene id
#             - type = gene
#             - size = 10
#             _Transcript_
#             - node id = ensembl transcript id
#             - type = transcript
#             - size = 5
#             - other attr: length of UTR sequence (for filtering?)
# 
#     NB: does each ensembl transcript id have only 1 utr? i.e. if the utr changes do you get a new transript id?

# <markdowncell>

# ### initiate graph object

# <codecell>

# undirected, no multi-edges
# latter point is important in adding edges to graph object
# in nx.Graph() node pair can only be connected once
# in nx.MultiGraph() they can have multiple edges
G = nx.Graph()

# collect node ids in dict for intitialisation and lookups
nodes = {
         # get unique set of genes and transcripts
         'genes': list(set([ utr[k]['gene'] for k in utr.keys()])), 
         'transcripts' : list(set([ utr[k]['transcript'] for k in utr.keys()]))
         }

# make sure lists are sorted - important for adding edges later
nodes['genes'].sort()
nodes['transcripts'].sort()

# check numbers
for k in nodes.keys():
    print(str(len(nodes[k])) + ' node ids collected for ' + k)

# <codecell>

# load gene nodes as
# node, attribute dict tuples list
G.add_nodes_from([
                  (k, {'type' : 'gene', 'size' : 200, '_id' : k, 'label' : k, 'col' : 'w'}) 
                  for k in nodes['genes']
                  ])

# <codecell>

# load UTR (transcript) nodes
G.add_nodes_from([
                  (k, {'type' : 'transcript', 'size' : 100, '_id' : k, 'label' : k, 'col' : 'b'}) 
                  for k in nodes['transcripts']
                  ]) 

# <codecell>

# add edges gene --> transcript
G.add_edges_from([(utr[k]['gene'], utr[k]['transcript'], {'type':'gene_transcript'}) for k in utr.keys()])
if VERBOSE:
    print(G.edges(data=True)[:3])

# <markdowncell>

# ### add edges to graph (with thresholding)

# <codecell>

# add edges from dataframe results
for row in results.iterrows():
    # ('to', {'from', 'weight', 'j_idx', 'j_dist'})
    
    if (row[1]['weight'] > 10 and row[1]['j_idx'] > 0.05):
        # add new edge with current weight
        G.add_edge(row[0], row[1]['from'], {'weight': row[1]['weight'], 'j_idx': row[1]['j_idx'], 'type': 'cocitation'})

print(str(len([i for i in G.edges(data=True) if i[2]['type']== 'cocitation' and i[2]['weight'] > 0])) + ' cocitation edges added')
if VERBOSE:
    print(G.edges(data=True)[:3])

# <codecell>

# add degree as node attribute
nx.set_node_attributes(G,'degree',{key:value for key, value in G.degree().items()})
nx.set_node_attributes(G,'weighted_degree',{key:value for key, value in G.degree(weight='weight').items()})

# <markdowncell>

# ## Plotting of results

# <codecell>

print('Analysing and plotting UTR-UTR cocitation graph.')

# <markdowncell>

# ### Genericised plotting function for graph objects

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

# <markdowncell>

# ### test plotting

# <codecell>

if VERBOSE and TEST:
    #generic_graph_plot(G)
    weighted_graph_plot(G, threshold_property='weight', threshold=100)
    weighted_graph_plot(G, threshold_property='j_idx', threshold=0.5)

# <markdowncell>

# ### plot graph for several thresholds

# <codecell>

if (TEST):
    for thresh in [i/100 for i in range(15,45, 5)]:
        p = weighted_graph_plot(G, threshold_property='j_idx', threshold=thresh)
        p.suptitle('Jaccard index threshold: ' + str(thresh), fontsize=25)
        p.savefig(OUTPUT_DIR + str(NMER_LEN) + 'mer_cocit.JI_gt_' + str(thresh) + '_.png')

# <markdowncell>

# ## save cocitation for later use

# <codecell>

nx.write_gml(G, OUTPUT_DIR + 'hsap.3utr.' + str(NMER_LEN) + 'mer_cocit.gml.gz') # automajickally compressed
# G = nx.read_gml(OUTPUT_DIR + 'hsap.3utr.cocit.gml', relabel=True)

# <codecell>


