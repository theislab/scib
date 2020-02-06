def Hbeta(D_row, beta):
    P = np.exp(- D_row * beta)
    sumP = P.sum()
    if (sumP == 0):
        H = 0    
        P = np.zeros(len(D_row))
    else:
        H = np.log(sumP) + beta * np.dot(D_row, P) / sumP
        P /= sumP
    return H, P

def compute_simpson_index(D = None, knn_idx = None, batch_labels = None, n_batches = None,
                          perplexity = 15, tol = 1e-5): 
    n = D.shape[0]
    P = np.zeros(D.shape[1])
    simpson = np.zeros(n)
    logU = np.log(perplexity)
    
    #loop over all cells
    for i in np.arange(0, n, 1):
        beta = 1
        # negative infinity
        betamin = -np.inf
        # positive infinity
        betamax = np.inf
        #get active row of D
        D_act = D[i,:]
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
                if (betamin== -np.inf):
                    beta /= 2
                else:
                    beta = (beta + betamin) / 2
    
      
            H, P = Hbeta(D_act, beta)
            Hdiff = H - logU
            tries += 1 
        
        if (H == 0):
            simpson[i] = -1
            continue
        
    
        #then compute Simpson's Index
        for b in np.arange(0, n_batches,1):
            q = np.flatnonzero(batch_labels[knn_idx[i]] == b) #indices of cells belonging to batch (b)
            if (len(q) > 0):
                sumP = np.sum(P[q])
                simpson[i] += sumP ** 2         
  
    return simpson
