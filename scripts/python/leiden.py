import igraph as ig
import leidenalg as la
import pandas as pd

protein_data = pd.read_csv("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.csv", index_col = [0])

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

sns.set(rc = {'figure.figsize':(10, 8)}, font_scale = 1)
sns.distplot(protein_data.mean(axis = 0), bins = 100)
plt.title('IP proteomics', fontsize = 25)
plt.xlabel('Mean Protein Exp', fontsize = 12)
plt.savefig("/home/degan/ip_proteomics/figures/UMAP/protein_mean.pdf")

import igraph as ig
from sklearn.neighbors import kneighbors_graph

## Building KNN graph
################################################################################

n_neighbor = 50
n_neighbor = 30
protein_data = protein_data.T
protein_ids = protein_data.index

knn_scRNAseq = kneighbors_graph(protein_data.values, n_neighbor, metric = 'euclidean', 
                                mode = 'connectivity').toarray()

knn_scRNAseq = pd.DataFrame(knn_scRNAseq, columns = protein_ids, index = protein_ids)
print(knn_scRNAseq.head())
print('Dimensions of knn: ' + str(knn_scRNAseq.shape))

knn_scRNAseq = knn_scRNAseq.stack().reset_index()
knn_scRNAseq = knn_scRNAseq.rename(columns = {'level_0': 'protein1', 'level_1': 'protein2', 
                                              0: 'connectivity'})

                                              
knn_scRNAseq = knn_scRNAseq.loc[knn_scRNAseq['connectivity'] != 0]

import cairocffi
from igraph import Graph, Plot
from IPython.display import Image
from igraph.drawing.text import TextDrawer

knn_scRNAseq = ig.Graph.TupleList([tuple(x) for x in knn_scRNAseq.values], 
                                  directed = False)
knn_scRNAseq.vs["label"] = knn_scRNAseq.vs['name']

visual_style = {}
visual_style["bbox"] = (800, 600)
visual_style["margin"] = 50

p = Plot("/home/degan/ip_proteomics/figures/UMAP/scRNAseq_graph.pdf", bbox = (800, 600), background = "white")
layout = knn_scRNAseq.layout_mds()

knn_scRNAseq.vs['color'] = ['red' if x['name']  in ['ESC_F04', 'ESC_F07'] 
                            else 'gold' for x in knn_scRNAseq.vs]
knn_scRNAseq.es['color'] = 'rgba(192, 192, 192, 0.3)'

p.add(knn_scRNAseq, layout = layout, vertex_size = 8, vertex_label_size = 1, 
      **visual_style)

p.redraw()
#ctx = cairocffi.Context(p.surface)
#ctx.set_font_size(20)
#drawer = TextDrawer(ctx, 'scRNAseq KNN Graph', halign = TextDrawer.CENTER)
#drawer.draw_at(300, 30, width = 200)
p.save()


import leidenalg
knn_consensus_cluster = leidenalg.find_partition(knn_scRNAseq, 
                                                 leidenalg.ModularityVertexPartition)
                                                 
                                                 
visual_style = {}
visual_style["bbox"] = (800, 600)
visual_style["margin"] = 50

p = Plot("/home/degan/ip_proteomics/figures/UMAP/consensus_graph.png", bbox = (800, 600), background = "white")
layout = knn_scRNAseq.layout_mds()

pal = ig.drawing.colors.ClusterColoringPalette(len(knn_consensus_cluster))
knn_scRNAseq.vs['color'] = pal.get_many(knn_consensus_cluster.membership)
p.add(knn_scRNAseq, layout = layout, vertex_size = 8, vertex_label_size = 1, 
      **visual_style)
p.redraw()

p.save()

################################################################################
## UMAP ####
################################################################################

from umap import UMAP
from sklearn.decomposition import PCA
from tqdm.auto import tqdm

protein_data = pd.read_csv("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.csv", index_col = [0])
protein_data = protein_data.T

X_scRNAseq = protein_data.values[:,0:(protein_data.shape[1]-1)]
Y_scRNAseq = protein_data.values[:,protein_data.shape[1]-1]

X_reduced = PCA(n_components = 20).fit_transform(X_scRNAseq)

model = UMAP(n_components = 2, min_dist = 0.3, n_neighbors = 60, 
             init = X_reduced[:, 0:2], n_epochs = 1000, verbose = 2)
umap = model.fit_transform(X_reduced)

plt.scatter(umap[:, 0], umap[:, 1], c = Y_scRNAseq, cmap = 'tab20', s = 2)
plt.title('UMAP: proteins', fontsize = 25); 
plt.xlabel("UMAP1", fontsize = 22); plt.ylabel("UMAP2", fontsize = 22)
plt.style.use('classic')
plt.savefig("/home/degan/ip_proteomics/figures/UMAP/protein_UMAP.pdf") 
plt.close()

                                                 
                                                 
                                            

