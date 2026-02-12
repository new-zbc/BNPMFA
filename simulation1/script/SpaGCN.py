import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt

folder = "simulation1/scenario1_7"
result = pd.DataFrame()

for data_index in range(10):
    data = str(data_index + 1)
    file_name = folder + "/"+ "SpaGCN/data/" + data +".h5ad"
    adata = sc.read_h5ad(file_name)
    adata.X = adata.X.astype("float32")

    x_pixel=adata.obs["col"].tolist()
    y_pixel=adata.obs["row"].tolist()
    #Calculate adjacent matrix
    s=1
    b=49

    
    #If histlogy image is not available, SpaGCN can calculate the adjacent matrix using the function below
    adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)

    #set hyper-parameters
    p=0.5 
    #Find the l value given p=
    l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)


    n_clusters = 3
    #print(n_clusters)
    #print(np.unique(obs_NA["Region"]))
    #Set seed
    r_seed=t_seed=n_seed=99
    #Seaech for suitable resolution
    res=spg.search_res(adata, adj, l, n_clusters, start=0.1, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)

    #run SpaGCN
    clf=spg.SpaGCN()
    clf.set_l(l)
    #Set seed
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    #Run
    clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
    y_pred, prob=clf.predict()
           

    adata.obs["pred"]= y_pred
    adata.obs["pred"]=adata.obs["pred"].astype('category')

    from sklearn import metrics
    obs_df = adata.obs.dropna()
    ari = metrics.adjusted_rand_score(obs_df['pred'], obs_df['label'])
    nmi = metrics.normalized_mutual_info_score(obs_df['pred'], obs_df['label'])
    ami = metrics.adjusted_mutual_info_score(obs_df['pred'], obs_df['label'])

    values= [ari, ami, nmi]
    result.append(values)

result_file = folder + "/summary/" + "SpaGCN.csv"
result.to_csv(result_file )
        