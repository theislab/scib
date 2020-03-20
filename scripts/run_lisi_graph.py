import warnings
warnings.filterwarnings('ignore')

import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import scipy.sparse as sparse
import networkx as nx
from networkx.algorithms.shortest_paths.weighted import all_pairs_dijkstra_path_length

import itertools
from multiprocessing import Pool
import multiprocessing

def consume(iterator, n=None):
    "Advance the iterator n-steps ahead. If n is None, consume entirely."
    # Use functions that consume iterators at C speed.
    if n is None:
        # feed the entire iterator into a zero-length deque
        collections.deque(iterator, maxlen=0)
    else:
        # advance to the empty slice starting at position n
        next(itertools.islice(iterator, n, n), None)

def take(n, iterable):
    "Return first n items of the iterable as a list"
    result = list(itertools.islice(iterable, n))
    return result

def Hbeta(D_row, beta):
    """
    Helper function for simpson index computation
    """
    P = np.exp(- D_row * beta)
    sumP = np.nansum(P)        
    if (sumP == 0):
            H = 0    
            P = np.zeros(len(D_row))
    else:
        H = np.log(sumP) + beta * np.nansum(D_row*P) / sumP
        P /= sumP
    return H, P

#LISI core function for shortest paths 
def compute_simpson_index_graph(D = None, batch_labels = None, n_batches = None, n_neighbors = 45,
                                  perplexity = 15, subsample = None, n_chunks = 10, chunk_no = 1,tol = 1e-5, verbose = False):
    """
    Simpson index of batch labels subsetted for each group.
    params:
        D: graph object
        batch_labels: a vector of length n_cells with batch info
        n_batches: number of unique batch labels
        n_neighbors: number of nearest neighbors
        perplexity: effective neighborhood size
        tol: a tolerance for testing effective neighborhood size
    returns:
        simpson: the simpson index for the neighborhood of each cell
    """
    #compute shortest paths of everyone to everyone
    dist = nx.all_pairs_dijkstra_path_length(D)
    n = len(batch_labels)
    P = np.zeros(n_neighbors)
    logU = np.log(perplexity)
    
    if subsample is None:
        subset = np.arange(0, n)
    else:
        subset = subsample

    #prepare chunk
    n_ch = n_chunks #number of chunks
    #get start and endpoint of chunk
    bounds = np.arange(0,n, np.ceil(n/n_ch).astype('int'))
    if chunk_no < n_ch - 1:
        chunk_ids = np.arange(bounds[chunk_no], bounds[chunk_no+1])
        if verbose:
            print(f"Entering chunk {chunk_no}.")
    else: #last chunk
        chunk_ids = np.arange(bounds[chunk_no], n)
        if verbose:
            print("Entering last chunk.")
                                            
    simpson = np.zeros(len(chunk_ids))
    #chunk has a start and an end
    if chunk_ids[0] != 0:
        consume(dist, chunk_ids[0]) #fast forward to first element of chunk
    
    #loop over all cells in chunk number
    for i in enumerate(chunk_ids): 
        if np.in1d(i[1], subset):
            #get neighbors and distances
            res = take(1,dist)
        else: #skip element
            consume(dist,1)
            continue
        #get sorted list of neighbours (keys) and distances (values)
        keys = np.array(list(res[0][1].keys()))
        values = np.array(list(res[0][1].values()))
        if len(keys)<n_neighbors:
            #not enough neighbors
            simpson[i[0]] = np.nan
            continue
                                                                                                                    
        beta = 1
        # negative infinity
        betamin = -np.inf
        # positive infinity
        betamax = np.inf
        #set distances
        D_act = values[1:][:n_neighbors]
        H, P = Hbeta(D_act, beta)
        Hdiff = H - logU
        tries = 0
        #first get neighbor probabilities
        while (np.logical_and(np.abs(Hdiff) > tol, tries < 50)):
            if (Hdiff > 0):
                betamin = beta
                if (betamax == np.inf): 
                    beta *= 2
                else:
                    beta = (beta + betamax) / 2
            else:
                betamax = beta
                if (betamin == -np.inf):
                    beta /= 2
                else:
                    beta = (beta + betamin) / 2
                
            H, P = Hbeta(D_act, beta)
            Hdiff = H - logU
            tries += 1
        
        if (H == 0):
            simpson[i[0]] = -1
            continue        
        
        #then compute Simpson's Index
        for b in np.arange(0, n_batches,1):
            knn_idx = keys[1:][:n_neighbors]
            non_nan_knn = knn_idx[np.invert(np.isnan(knn_idx))].astype('int')
            q = np.flatnonzero(batch_labels[non_nan_knn] == b) #indices of cells belonging to batch (b)
            if (len(q) > 0):
                sumP = np.sum(P[q])
                simpson[i[0]] += sumP ** 2         
    
    #return only entries in the subset
    simpson = simpson[np.in1d(chunk_ids,subset)] 

    return simpson

def lisi_graph_py(adata, batch_key, n_neighbors = 90, perplexity=None, subsample = None, chunk_id = None, 
                  n_chunks = 10, verbose=False):
    """
    Compute LISI score on shortes path based on kNN graph provided in the adata object. 
    By default, perplexity is chosen as 1/3 * number of nearest neighbours in the knn-graph.
    """
    
    if 'neighbors' not in adata.uns:
        raise AttributeError(f"key 'neighbors' not found. Please make sure that a " +
                              "kNN graph has been computed")    
    elif verbose:                                                    
        print("using precomputed kNN graph")
                                                                        
    #get knn index matrix
    if verbose:
        print("Convert nearest neighbor matrix and distances for LISI.")
               
    batch = adata.obs[batch_key].cat.codes.values        
    n_batches = len(np.unique(adata.obs[batch_key])) 
                                        
    if perplexity is None or perplexity >=n_neighbors:
        # use LISI default
        perplexity = np.floor(n_neighbors/3)
                                                                                                                                
    # run LISI in python
    if verbose:
        print("Compute shortest paths") 
                                                                                                
    #turn connectivities matrix into graph
    G = nx.from_scipy_sparse_matrix(adata.uns['neighbors']['connectivities'])  
    
    if verbose:
        print("LISI score estimation")
    
    if chunk_id is None:
        simpson_estimate_dict = {}
        for idx in np.arange(0,n_chunks):
            tmp_res = compute_simpson_index_graph(D = G, 
                                                  batch_labels = batch,                           
                                                  n_batches = n_batches,
                                                  perplexity = perplexity, 
                                                  subsample = subsample,
                                                  n_neighbors = n_neighbors,
                                                  n_chunks = n_chunks,
                                                  chunk_no = idx,
                                                  verbose = True
                                                 )
        #create array from dictionary
        simpson_estimate_batch = np.concatenate([simpson_estimate_dict[ids] for ids in simpson_estimate_dict.keys()])
    else:
        simpson_estimate_batch = compute_simpson_index_graph(D = G, 
                                                             batch_labels = batch,
                                                             n_batches = n_batches,
                                                             perplexity = perplexity,
                                                             n_neighbors = n_neighbors,
                                                             subsample = subsample,
                                                             n_chunks = n_chunks,
                                                             chunk_no = chunk_id
                                                            )
    simpson_est_batch = 1/simpson_estimate_batch
    # extract results
    d = {batch_key : simpson_est_batch}
    lisi_estimate = pd.DataFrame(data=d, index=np.arange(0,len(simpson_est_batch)))
    
    return lisi_estimate

# types of integration output
RESULT_TYPES = [
    "full", # reconstructed expression data
    "embed", # embedded/latent space
    "knn" # only corrected neighbourhood graph as output
]


if __name__=='__main__':
    """
    read adata object, compute shortest paths using the dijkstra algorithm and compute iLISI score.
    """

    import argparse
    import os

    parser = argparse.ArgumentParser(description='Compute LISI graph distances based on shortest paths.')

    parser.add_argument('-i', '--input', required = True)
    #parser.add_argument('-o', '--output', required = True, help = 'output directory')
    parser.add_argument('-v', '--verbose')
    parser.add_argument('-t', '--type', required=True, choices=RESULT_TYPES, 
            help='Type of result: full, embed, knn\n full: scanorama, seurat, MNN\n embed: scanorama, Harmony\n knn: BBKNN')
    parser.add_argument('-b', '--batch_key', required = True)
    parser.add_argument('-n', '--nodes', default = None, help='Number of nodes for multiprocessing.')
    parser.add_argument('-s', '--subsample', type=float, default = None, help='Fraction to which should be subsampled.')

    args = parser.parse_args()

    verbose = args.verbose
    type_ = args.type
    batch_key = args.batch_key
    nodes = args.nodes
    subsample = args.subsample
    base = os.path.basename(args.input).split('.h5ad')[0]
    outdir = os.path.split(args.input)[0]

    print("reading adata input file")
    if os.stat(args.input).st_size>0:
        adata = sc.read(args.input, cache = True)
        sc.pp.subsample(adata, n_obs=2000)
        if (type_ == 'embed'):
            sc.pp.neighbors(adata, use_rep = 'X_emb')
        if (type_ == 'full'):
            sc.pp.pca(adata, svd_solver = 'arpack')
            sc.pp.neighbors(adata)
        #if knn - do not compute a new neighbourhood graph (it exists already)
        
        #prepare simpson call
        batch = adata.obs[batch_key].cat.codes.values
        n_batches = len(np.unique(adata.obs[batch_key]))
        n_neighbors = 90
        perplexity = np.floor(n_neighbors/3)
        G = nx.from_scipy_sparse_matrix(adata.uns['neighbors']['connectivities'])
        
        if subsample is not None:
            subset = np.random.choice(np.arange(0,adata.n_obs), 
                     np.floor(subsample*adata.n_obs).astype('int'),
                     replace=False
                     )
        else:
            subset = None


        #set up multiprocessing
        if nodes is None:
            n_processes = np.max([multiprocessing.cpu_count() - 5, 
                               np.ceil(multiprocessing.cpu_count()/2)]).astype('int')
        else:
            n_processes = nodes

        print(f"{n_processes} processes started.")
        pool = Pool(processes=n_processes)
        count = np.arange(0, n_processes)
        
        #create argument list for each worker
        results = pool.starmap(compute_simpson_index_graph, zip(itertools.repeat(G),
                                                                itertools.repeat(batch),
                                                                itertools.repeat(n_batches),
                                                                itertools.repeat(n_neighbors),
                                                                itertools.repeat(perplexity),
                                                                itertools.repeat(subset),
                                                                itertools.repeat(n_processes),
                                                                count))
        pool.close()
        pool.join()
        #q = multiprocessing.Queue()
        #processes = []
        #results = []
        #for i in count:
        #    p = multiprocessing.Process(target=lisi_graph_py, args=(adata, batch_key, 90, None, i, n_processes))
        #    processes.append(p)
        #    p.start()
        #for p in processes:
        #    result = q.get()
        #    results.append(results)
        #for p in processes:
        #    p.join()
        #lisi_result = lisi_graph_py(adata=adata, batch_key = batch_key, n_neighbors = 90)
        #print(results)
        simpson_est_batch = 1/np.concatenate(results)
        #print(np.nanmedian(simpson_est_batch))
        # extract results
        d = {batch_key : simpson_est_batch}
        if subset is None:
            simpson_index = pd.DataFrame(data=d, index=np.arange(0,len(simpson_est_batch)))
        else:
            simpson_index = pd.DataFrame(data=d, index=np.sort(subset))
        simpson_index.to_csv(os.path.join(outdir, f'{base}_{type_}_lisi_graph.csv'))
    else:
        print("No file found. Doing nothing.")

