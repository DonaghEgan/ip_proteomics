################################################################################
## UMAP ####
################################################################################

from igraph import Graph, Plot
from IPython.display import Image
from igraph.drawing.text import TextDrawer
import leidenalg
from umap import UMAP
from sklearn.decomposition import PCA
from tqdm.auto import tqdm
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import igraph as ig
import pandas as pd

protein_data = pd.read_csv("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.csv", index_col = [0])
protein_data = protein_data.T

X_scRNAseq = protein_data.values[:,0:(protein_data.shape[1]-1)]
Y_scRNAseq = protein_data.values[:,protein_data.shape[1]-1]

X_reduced = PCA(n_components = 4).fit_transform(X_scRNAseq)

model = UMAP(n_components = 4, min_dist = 0.3, n_neighbors = 70, 
             init = X_reduced[:, 0:4], n_epochs = 1000, verbose = 2)
umap = model.fit_transform(X_reduced)

plt.scatter(umap[:, 0], umap[:, 1], c = Y_scRNAseq, cmap = 'tab20', s = 8)
plt.title('UMAP: IP', fontsize = 16); 
plt.xlabel("UMAP1", fontsize = 12); plt.ylabel("UMAP2", fontsize = 12)
plt.style.use('classic')
plt.savefig("/home/degan/ip_proteomics/figures/UMAP/protein_UMAP.pdf") 
plt.close()
