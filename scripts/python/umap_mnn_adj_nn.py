## modules ####
################################################################################

import pandas as pd
import numpy as np
import pandas as pd
import networkx as nx
import scipy
import numba
import pickle
import gzip
from umap.umap_ import nearest_neighbors
from matplotlib import cm
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt


## Reading in data ####
################################################################################

X = pd.read_csv("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.csv", index_col = [0])
X = X.T
Y_scRNAseq = X.values[:,X.shape[1]-1]

## Parameters ####
################################################################################

# Random State
from sklearn.utils import check_random_state
random_state = check_random_state(0)

# Distance Metric to Use
metric = 'euclidean'

# number of neighbors for computing k-neighbor graph
n_neighbors = 30

# new number of neighbors to search for
new_n_neighbors = 30

## Calculate a,b hyperparams given the min_dist ####
################################################################################

from umap.umap_ import find_ab_params

min_dist = 0.1
a, b = find_ab_params(1.0, min_dist)

## Calculate the weighted nearest neighbor graph ####
################################################################################

knn_indices, knn_dists, knn_search_index = nearest_neighbors(
    X,
    n_neighbors=n_neighbors,
    metric = metric,
    metric_kwds = {},
    angular=False,
    random_state = random_state,
    low_memory=True,
    use_pynndescent=True,
    n_jobs=1,
    verbose=True,
)

## Calculate min spanning tree ####
################################################################################

def min_spanning_tree(knn_indices, knn_dists, n_neighbors, threshold):
  
  rows = np.zeros(knn_indices.shape[0] * n_neighbors, dtype=np.int32)
  cols = np.zeros(knn_indices.shape[0] * n_neighbors, dtype=np.int32)
  vals = np.zeros(knn_indices.shape[0] * n_neighbors, dtype=np.float32)
  
  pos = 0
  for i, indices in enumerate(knn_indices):
    for j, index in enumerate(indices[:threshold]):
      if index == -1:
        continue
      rows[pos] = i 
      cols[pos] = index
      vals[pos] = knn_dists[i][j]
      pos += 1
  
  matrix = scipy.sparse.csr_matrix((vals, (rows, cols)), shape=(knn_indices.shape[0], knn_indices.shape[0]))
  Tcsr = scipy.sparse.csgraph.minimum_spanning_tree(matrix)
  
  Tcsr = scipy.sparse.coo_matrix(Tcsr)
  weights_tuples = zip(Tcsr.row, Tcsr.col, Tcsr.data)
  

  sorted_weights_tuples = sorted(weights_tuples, key=lambda tup: tup[2])
  
  return sorted_weights_tuples 

## Connected mnn
################################################################################

import copy
import heapq

def create_connected_graph(mutual_nn, total_mutual_nn, knn_indices, knn_dists, n_neighbors, connectivity):
  connected_mnn = copy.deepcopy(mutual_nn)
  
  if connectivity == "nearest":
    for i in range(len(knn_indices)): 
      if len(mutual_nn[i]) == 0:
        first_nn = knn_indices[i][1]
        if first_nn != -1:
          connected_mnn[i].add(first_nn) 
          connected_mnn[first_nn].add(i) 
          total_mutual_nn += 1
    return connected_mnn

            
      
  #Create graph for mutual NN
  rows = np.zeros(total_mutual_nn, dtype=np.int32)
  cols = np.zeros(total_mutual_nn, dtype=np.int32)
  vals = np.zeros(total_mutual_nn, dtype=np.float32)
  pos = 0
  for i in connected_mnn:
    for j in connected_mnn[i]:
      rows[pos] = i 
      cols[pos] = j
      vals[pos] = 1
      pos += 1
  graph = scipy.sparse.csr_matrix((vals, (rows, cols)), shape=(knn_indices.shape[0], knn_indices.shape[0]))

  #Find number of connected components
  n_components, labels = scipy.sparse.csgraph.connected_components(csgraph=graph, directed=True, return_labels=True, connection= 'strong')
  print(n_components)
  label_mapping = {i:[] for i in range(n_components)}

  for index, component in enumerate(labels):
    label_mapping[component].append(index)



  #Find the min spanning tree with KNN
  sorted_weights_tuples = min_spanning_tree(knn_indices, knn_dists, n_neighbors, n_neighbors)
  

  #Add edges until graph is connected
  for pos,(i,j,v) in enumerate(sorted_weights_tuples):

    if connectivity == "full_tree":
      connected_mnn[i].add(j)
      connected_mnn[j].add(i) 
      
      
    elif connectivity == "min_tree" and labels[i] != labels[j]:
      if len(label_mapping[labels[i]]) < len(label_mapping[labels[j]]):
        i, j = j, i
        
      connected_mnn[i].add(j)
      connected_mnn[j].add(i)
      j_pos = label_mapping[labels[j]]
      labels[j_pos] = labels[i]
      label_mapping[labels[i]].extend(j_pos)

  return connected_mnn  

##
################################################################################

#Search to find adjacent neighbors 
def find_new_nn(knn_indices, knn_dists, knn_indices_pos, connected_mnn, n_neighbors_max, verbose=False):
  
  new_knn_dists= [] 
  new_knn_indices = []
  
  for i in range(len(knn_indices)): 
    #print(i)
    min_distances = []
    min_indices = []
    #Initialize vars
    heap = [(0,0,i)]
    mapping = {}
          
    seen = set()
    heapq.heapify(heap) 
    while(len(min_distances) < n_neighbors_max and len(heap) >0):
      dist, hop, nn = heapq.heappop(heap)
      if nn == -1:
        continue
      #For adjacent, only considering one hop away
      if nn not in seen and hop <= 1:
        min_distances.append(dist)
        min_indices.append(nn)
        seen.add(nn)
        neighbor = connected_mnn[nn]
        
        for nn_nn in neighbor:
          if nn_nn not in seen and hop <= 0:
            distance = 0
            if nn_nn in knn_indices_pos[nn]:
              pos = knn_indices_pos[nn][nn_nn]
              distance = knn_dists[nn][pos] 
            else:
              pos = knn_indices_pos[nn_nn][nn]
              distance = knn_dists[nn_nn][pos] 
            distance += dist
            
            if nn_nn not in mapping:
              mapping[nn_nn] = distance
              heapq.heappush(heap, (distance, hop+1, nn_nn))
            elif mapping[nn_nn] > distance:
              mapping[nn_nn] = distance
              heapq.heappush(heap, (distance, hop+1, nn_nn))
    
    if len(min_distances) < n_neighbors_max:
      for j in range(n_neighbors_max-len(min_distances)):
        min_indices.append(-1)
        min_distances.append(np.inf)
    
    new_knn_dists.append(min_distances)
    new_knn_indices.append(min_indices)
    
    if verbose and i % int(len(knn_dists) / 10) == 0:
      print("\tcompleted ", i, " / ", len(knn_dists), "epochs")
  return new_knn_dists, new_knn_indices



#Calculate the connected mutual nn graph
def mutual_nn_nearest(knn_indices, knn_dists, n_neighbors, n_neighbors_max, connectivity="min_tree", verbose=False):
  mutual_nn = {}
  nearest_n= {}

  knn_indices_pos = [None] * len(knn_indices)

  total = 0
  
  for i, top_vals in enumerate(knn_indices):
    nearest_n[i] = set(top_vals)
    knn_indices_pos[i] = {}
    for pos, nn in enumerate(top_vals):
      knn_indices_pos[i][nn] = pos
  
  total_mutual_nn = 0
  for i, top_vals in enumerate(knn_indices):
    mutual_nn[i] = set()
    for ind, nn in enumerate(top_vals):
      if nn != -1 and (i in nearest_n[nn] and i != nn):
        mutual_nn[i].add(nn)
        total_mutual_nn += 1

  
  connected_mnn = create_connected_graph(mutual_nn, total_mutual_nn, knn_indices, knn_dists, n_neighbors, connectivity )
  new_knn_dists, new_knn_indices = find_new_nn(knn_indices, knn_dists, knn_indices_pos, connected_mnn, n_neighbors_max, verbose)

  
  return connected_mnn, mutual_nn, np.array(new_knn_indices), np.array(new_knn_dists)  

    
connected_mnn,mutual_nn, new_knn_indices, new_knn_dists  = mutual_nn_nearest(knn_indices, knn_dists, n_neighbors, 
new_n_neighbors, connectivity= "nearest" , verbose = True)

## Calculate the Fuzzy Simplicial Set
################################################################################

from umap.umap_ import compute_membership_strengths

INT32_MIN = np.iinfo(np.int32).min + 1
INT32_MAX = np.iinfo(np.int32).max - 1

SMOOTH_K_TOLERANCE = 1e-5
MIN_K_DIST_SCALE = 1e-3
NPY_INFINITY = np.inf

@numba.njit(
    locals={
        "psum": numba.types.float32,
        "lo": numba.types.float32,
        "mid": numba.types.float32,
        "hi": numba.types.float32,
    },
    fastmath=True,
)  

#Modified to use a variable k as explained paper
def smooth_knn_dist(distances, k, n_iter=64, local_connectivity=1.0, bandwidth=1.0):

    #target = np.log2(k) * bandwidth
    rho = np.zeros(distances.shape[0], dtype=np.float32)
    result = np.zeros(distances.shape[0], dtype=np.float32)

    mean_distances = np.mean(distances)

    for i in range(distances.shape[0]):

        
        lo = 0.0
        hi = NPY_INFINITY
        mid = 1.0

        # TODO: This is very inefficient, but will do for now. FIXME
        ith_distances = distances[i]
        non_zero_dists = ith_distances[ith_distances > 0.0]
        
        # New k
        non_inf_dists = ith_distances[ith_distances != np.inf]
        target = np.log2(len(non_inf_dists)) * bandwidth
        #print(target)
        
        
        if non_zero_dists.shape[0] >= local_connectivity:
            index = int(np.floor(local_connectivity))
            interpolation = local_connectivity - index
            if index > 0:
                rho[i] = non_zero_dists[index - 1]
                if interpolation > SMOOTH_K_TOLERANCE:
                    rho[i] += interpolation * (
                        non_zero_dists[index] - non_zero_dists[index - 1]
                    )
            else:
                rho[i] = interpolation * non_zero_dists[0]
        elif non_zero_dists.shape[0] > 0:
            rho[i] = np.max(non_zero_dists)

        for n in range(n_iter):

            psum = 0.0
            for j in range(1, distances.shape[1]):
                d = distances[i, j] - rho[i]
                if d > 0:
                    psum += np.exp(-(d / mid))
                else:
                    psum += 1.0

            if np.fabs(psum - target) < SMOOTH_K_TOLERANCE:
                break

            if psum > target:
                hi = mid
                mid = (lo + hi) / 2.0
            else:
                lo = mid
                if hi == NPY_INFINITY:
                    mid *= 2
                else:
                    mid = (lo + hi) / 2.0

        result[i] = mid

        # TODO: This is very inefficient, but will do for now. FIXME
        if rho[i] > 0.0:
            mean_ith_distances = np.mean(ith_distances)
            if result[i] < MIN_K_DIST_SCALE * mean_ith_distances:
                result[i] = MIN_K_DIST_SCALE * mean_ith_distances
        else:
            if result[i] < MIN_K_DIST_SCALE * mean_distances:
                result[i] = MIN_K_DIST_SCALE * mean_distances

    return result, rho

#Calculate fuzzy_simplicial_set with new variable k
def fuzzy_simplicial_set(
    X,
    n_neighbors,
    random_state,
    metric,
    metric_kwds={},
    knn_indices=None,
    knn_dists=None,
    angular=False,
    set_op_mix_ratio=1.0,
    local_connectivity=1.0,
    apply_set_operations=True,
    verbose=False,
    return_dists=None,
):
    if knn_indices is None or knn_dists is None:
        knn_indices, knn_dists, _ = nearest_neighbors(
            X,
            n_neighbors,
            metric,
            metric_kwds,
            angular,
            random_state,
            verbose=verbose,
        )

    knn_dists = knn_dists.astype(np.float32)

    sigmas, rhos = smooth_knn_dist(
        knn_dists,
        float(n_neighbors),
        local_connectivity=float(local_connectivity),
    )

    rows, cols, vals, dists = compute_membership_strengths(
        knn_indices, knn_dists, sigmas, rhos, return_dists
    )

    result = scipy.sparse.coo_matrix(
        (vals, (rows, cols)), shape=(X.shape[0], X.shape[0])
    )
    result.eliminate_zeros()

    if apply_set_operations:
        transpose = result.transpose()

        prod_matrix = result.multiply(transpose)

        result = (
            set_op_mix_ratio * (result + transpose - prod_matrix)
            + (1.0 - set_op_mix_ratio) * prod_matrix
        )

    result.eliminate_zeros()

    if return_dists is None:
        return result, sigmas, rhos
    else:
        if return_dists:
            dmat = scipy.sparse.coo_matrix(
                (dists, (rows, cols)), shape=(X.shape[0], X.shape[0])
            )

            dists = dmat.maximum(dmat.transpose()).todok()
        else:
            dists = None

        return result, sigmas, rhos, dists  
      
# build fuzzy_simplicial_set
P, sigmas, rhos = fuzzy_simplicial_set(
    X = X,
    n_neighbors = new_n_neighbors,
    metric = metric,
    random_state = random_state,
    knn_indices= new_knn_indices,
    knn_dists = new_knn_dists,
)


## Find the low dimensional representation
#################################################################################

from umap.umap_ import simplicial_set_embedding
from umap.umap_ import dist
#Dimensionality of the low dimensional representation
n_components = 2

negative_sample_rate = 5

embeddings , aux_data = simplicial_set_embedding(
    data = X,
    graph = P,
    n_components = n_components,
    initial_alpha = 1.0,
    a = a,
    b = b,
    gamma = 1.0,
    negative_sample_rate = negative_sample_rate,
    n_epochs = -1,
    init = "spectral",
    random_state = check_random_state(0),
    metric = metric,
    metric_kwds = {},
    densmap = False,
    densmap_kwds = {},
    output_dens = False,
    output_metric= dist.named_distances_with_gradients["euclidean"],
    output_metric_kwds={},
    euclidean_output=True,
    parallel=False,
    verbose=True,
)

#plot for 2d
colors = []
colors += cm.get_cmap("Set3").colors
colors += cm.get_cmap("Set2").colors
my_cmap = ListedColormap(colors)
#print(my_cmap)
plt.rcParams.update(plt.rcParamsDefault)
plt.scatter(embeddings[:, 0], embeddings[:, 1], c = Y_scRNAseq, cmap="Spectral",  s = 4)
plt.savefig("/home/degan/ip_proteomics/figures/UMAP/protein_UMAP_mnn.pdf") 
plt.close()


















