from math import log, exp
import numpy as np 
import scipy.sparse as sparse
import networkx as nx
cimport cython
DTYPE = np.float64


cdef extern from "math.h":
    float INFINITY

@cython.boundscheck(False)
@cython.wraparound(False)
def compute_simpson_index_graph_cy(double[::1] data, int[::1] indices, int[::1] indptr, 
                                long[::1] batch_labels, int n_batches, int n_neighbors,
                                int perplexity, long[::1] subsample, 
                                int n_chunks, int chunk_no, float tol 
                                ):
    """
    Simpson index of batch labels subsetted for each group.
    params:
        D: graph object (as sparse matrix)
        batch_labels: a vector of length n_cells with batch info
        n_batches: number of unique batch labels
        n_neighbors: number of nearest neighbors
        perplexity: effective neighborhood size
        tol: a tolerance for testing effective neighborhood size
    returns:
        simpson: the simpson index for the neighborhood of each cell
    """
    #compute shortest paths of everyone to everyone
    #Update: We don't actually need that, because we can compute 
    #the distance from one to all others when we actually need it
    #dist = nx.all_pairs_dijkstra_path_length(D)
    connectivity = sparse.csr_matrix((data, indices, indptr))
    D = nx.from_scipy_sparse_matrix(connectivity)
    
    
    cdef int n = len(batch_labels)

    cdef double[::1] P = np.zeros(n_neighbors, dtype=DTYPE)
    cdef double[::1] P_tmp = np.zeros(n_neighbors, dtype=DTYPE)
    cdef double logU = log(perplexity)
    
    cdef int subset_out = 1
    cdef int sub_len = len(subsample)
    cdef Py_ssize_t x
    cdef Py_ssize_t y
    cdef Py_ssize_t i 
    cdef Py_ssize_t q 
    cdef long c_id
    cdef long s_id
    cdef int k = 0
    cdef int j = 0
    

    cdef int chunk_len = n
    cdef int max_range = n
    cdef int min_range = 0
    cdef long[::1] chunk_ids = np.zeros(max_range-min_range, dtype=np.long)
    cdef dict res
    cdef long[::1] keys = np.zeros(n_neighbors, dtype=np.long)
    cdef double[::1] values = np.zeros(n_neighbors, dtype=DTYPE)
    cdef long key 
    cdef double value
    cdef double beta = 1
    #define infinity
    inf = INFINITY #if INFINITY != 0 else float('inf')
    
    cdef double betamin = -inf
    cdef double betamax = inf
    cdef double[::1] D_act = np.zeros(n_neighbors, dtype=DTYPE)
    cdef double H = 0
    cdef double Hdiff = 0 
    cdef double sumP_tmp = 0
    cdef double tmp_beta = 0
    cdef int tries = 0
    cdef int idx = 0
    cdef long[::1] knn_idx = np.zeros(n_neighbors, dtype=np.long)
    
    cdef long[::1] batch = np.zeros(n_neighbors, dtype=np.long)
    cdef double[:,::1] B = np.zeros((n_neighbors, n_batches), dtype=DTYPE) 
    cdef double sumP = 0
    cdef double p_i = 0
    cdef double B_i = 0 
    cdef double dot_sum = 0
    
    #prepare chunk
    if n_chunks>1:
        #get start and endpoint of chunk
        chunk_len = n//n_chunks
        max_range =(1+chunk_no)* chunk_len if chunk_no < n_chunks - 1 else n
        min_range = chunk_no*chunk_len
    
    if sub_len>0:
        subset_out = 0
        for y in range(sub_len):
            s_id = subsample[y]
            subset_out += (s_id >= min_range and s_id<max_range)
    
    cdef long[::1] tmp = np.zeros(subset_out, dtype=np.long)
    #remove chunk_ids, which are not in subsample
    if sub_len>0:
        for x in range(min_range, max_range):
            k = 0
            for y in range(sub_len):
                s_id = subsample[y]
                k+= (x==s_id)
            if (k>0):    
                tmp[j] = x
                j+=1
        chunk_ids = tmp
        chunk_len = subset_out
    else:
        for x in range(chunk_len):
            chunk_ids[x] = x + min_range
    
    simpson = np.zeros(chunk_len, dtype=DTYPE)
    cdef double[::1] result_view = simpson
    
    
    #loop over all cells in chunk number
    for i in range(chunk_len):
        #get neighbors and distances
        res = nx.single_source_dijkstra_path_length(D, chunk_ids[i])
        if len(res)<n_neighbors:
            #not enough neighbors
            result_view[i] = 1 # np.nan #set nan for testing
            continue
        #get sorted list of neighbours (keys) and distances (values)
        k = 0 
        for key, value in res.items():
            keys[k] = key
            values[k] = value
            k +=1
        #keys = np.array(list(res.keys()))
        #values = np.array(list(res.values()))
        
        #start lisi estimation
        beta = 1
        # negative infinity
        betamin = -inf
        # positive infinity
        betamax = inf
        #set distances
        D_act = values[1:][:n_neighbors]
        sumP_tmp = 0
        dot_sum = 0
        #replace Hbeta
        for x in range(n_neighbors):
            P_tmp[x] = exp(-D_act[x] * beta)
            sumP_tmp += P_tmp[x]
        if (sumP_tmp == 0):
            #H = 0   
            for x in range(n_neighbors):
                P_tmp[x] = 0
        else:
            for x in range(n_neighbors):
                tmp_beta = D_act[x]*P_tmp[x]
                dot_sum += tmp_beta
        
            H = log(sumP_tmp) + beta * dot_sum / sumP_tmp
            for x in range(n_neighbors):
                P[x] = P_tmp[x]/sumP_tmp
        #H, P = Hbeta(D_act, beta)
        Hdiff = H - logU
        tries = 0
        #first get neighbor probabilities
        while abs(Hdiff) > tol and tries < 50:
            if (Hdiff > 0):
                betamin = beta
                if (betamax == inf): 
                    beta *= 2
                else:
                    beta = (beta + betamax) / 2
            else:
                betamax = beta
                if (betamin == -inf):
                    beta /= 2
                else:
                    beta = (beta + betamin) / 2
                
            #replace H, P = Hbeta(D_act, beta)
            sumP_tmp = 0
            dot_sum = 0
            #replace Hbeta
            for x in range(n_neighbors):
                P_tmp[x] = exp(-D_act[x] * beta)
                sumP_tmp += P_tmp[x]
            if (sumP_tmp == 0):
                #H = 0   
                for x in range(n_neighbors):
                    P_tmp[x] = 0
            else:
                for x in range(n_neighbors):
                    tmp_beta = D_act[x]*P_tmp[x]
                    dot_sum += tmp_beta
        
                H = log(sumP_tmp) + beta * dot_sum / sumP_tmp
                for x in range(n_neighbors):
                    P[x] = P_tmp[x]/sumP_tmp
            #H, P = Hbeta(D_act, beta)        
            Hdiff = H - logU
            tries += 1
        
        if (H == 0):
            result_view[i] = -1
            continue 
        #then compute Simpson's Index
        knn_idx = keys[1:][:n_neighbors]
        for x in range(n_neighbors):
            batch[x] = batch_labels[knn_idx[x]] 
            for p in range(n_batches):
                B[x,p] = 0
            B[x, batch[x]] = 1    
        
        for q in range(n_batches):
            sumP = 0
            #replace matrix multiplication
            for x in range(n_neighbors):
                p_i = P[x]
                B_i = B[x,q]
                dot_sum = p_i * B_i
                sumP+= dot_sum
                
            result_view[i] += sumP ** 2 #sum squares
        
    return simpson 
