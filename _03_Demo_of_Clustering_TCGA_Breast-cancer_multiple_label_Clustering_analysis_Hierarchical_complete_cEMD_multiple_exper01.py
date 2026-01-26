####+++++++++++
# Notes: This is a demo python scripts of how we perform the Clustering, and calculated the corresponding metrics considering different etiology labels, 
# based on the prepared data of signature exposures in TCGA cancer samples.

# In this demo scripts, we performed the 'Hierarchical' clustering (with linkage = complete) for Breast-cancer samples based normalized exposures, 
# using the cEMD distance matrix, and the clustering is performed based on subset of samples (50%) for multiple time.

# We also varied cancer types, distance types, clustering algorithms for other custering cases, also based on sub-samples dataset for multiple times, in order verify the robustness of results
# Similar analysis also performed on Hartwig data 
####+++++++++++

# ++++++++++++++++
# import modules 
# ++++++++++++++++

from fileinput import filename
import os
import time
import random
import pandas as pd
import numpy as np
import pickle

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import v_measure_score
from sklearn.metrics import homogeneity_score, completeness_score

os.chdir('/path/to/EMD_analysis')
os.getcwd()

clusteringPath = '/path/to/EMD_analysis/..'

### Define the function for Hierarchical_complete Clustering
def hierarchical_clustering(distMat, n_clusters, linkage_method='complete'):
    '''
    Perform hierarchical clustering and return cluster labels
    
    distMat: The given distance matrix (must be in a format suitable for AgglomerativeClustering)
    n_clusters: The number of clusters to find
    linkage_method: The linkage criterion to use. Options include 'ward', 'complete', 'average', and 'single'
    
    Returns:
    cluster_labels: The labels of each point in the dataset after clustering
    '''
    # Create a model for clustering with the given number of clusters and linkage method
    # 'precomputed' metric is used because we are providing a precomputed distance matrix
    model = AgglomerativeClustering(n_clusters=n_clusters, 
                                    metric='precomputed',
                                    linkage=linkage_method)

    # Perform the clustering and extract the cluster labels
    cluster_labels = model.fit_predict(distMat)

    # Return the cluster labels
    return cluster_labels   
    
### Load the related data
cancer_types = ['Breast-cancer']
distTypes = ['cEMD'] 
labelTypes = ['label_Apobec', 'label_MMR', 'label_Tobacco', 'label_UV', 'label_POLE', 'label_ClockLike', 'label_BER', 'label_Platinum'] 

# define the dict to hold the scores, using cancer type as keys
dict_allScores_hierarchical = {}

for cancer_type in cancer_types:
    
    ## Load the data
    with open(f'TCGA_all_distDFs_dict_for_{cancer_type}_Sample_exper01.pickle', 'rb') as handle:
        dict_distDF = pickle.load(handle)
    df_frac = pd.read_csv(f'_TCGA {cancer_type} Mutational Signatures Fraction Multi-Label Data_exper01.csv', index_col=0)

    start_time = time.time() # The start time
    
    ## Perform the random subsampling, and perform clustering on sub-sampled samples
    for round in range(10):
        
        # Randomly sample 50% of the samples based on the index
        sampled_df = df_frac.sample(frac=0.5, random_state=round)
        
        # index of sub-samples
        index = sampled_df.index
        
        # Define the max number of clusters
        max_clusters = len(sampled_df)
    
        ## Dicts for holding output of different scores, the keys will be each distance metrices, 
        ## and the values will be the dataframe of scores
        dict_df_vScores = dict() ## dict of data frames for V scores,
        dict_df_homogeneityScores = dict() ## dict of data frames for Homogeneity, 
        dict_df_completenessScores = dict() ## dict of data frames for Completeness, 
    
        ## Perform the computation
        for distType in distTypes:
    
            # Taking the subset of distance matrix
            distMat = dict_distDF[distType]
            sub_distMat = distMat.loc[index, index]
            
            # Define the dataframes to accept results
            vScore_df = pd.DataFrame(index=range(2, max_clusters, 1), columns=labelTypes)
            homogeneityScore_df = pd.DataFrame(index=range(2, max_clusters, 1), columns=labelTypes)
            completenessScore_df = pd.DataFrame(index=range(2, max_clusters, 1), columns=labelTypes)
    
            # Looping the number of clusters
            clustering_columns = []
            for n_cluster in range(2, max_clusters, 1):
                
                cluster_label = hierarchical_clustering(sub_distMat, n_cluster)
                cluster_name = 'Cluster_K='+str(n_cluster)
                clustering_columns.append(pd.Series(cluster_label, name=cluster_name, index=sampled_df.index))
                
                # Looping the ground label types: 
                for labelType in labelTypes:
    
                    ## Taking the label type
                    ground_label = sampled_df[labelType]
                    
                    ## V-Score
                    v_score = v_measure_score(ground_label, cluster_label)
                    ## Homogenetity
                    homogeneity = homogeneity_score(ground_label, cluster_label)
                    ## Completeness
                    completeness = completeness_score(ground_label, cluster_label)
    
                    # accept the results
                    vScore_df.loc[n_cluster, labelType] = v_score
                    homogeneityScore_df.loc[n_cluster, labelType] = homogeneity
                    completenessScore_df.loc[n_cluster, labelType] = completeness
            
            ## Save the clustering results
            filename = f'_TCGA {cancer_type} Mutational Signatures Fraction Multi-Label Data_{distType}_Hierarchical_complete_clustering_results_round{round}_exper01.csv'
            clusteringfile = os.path.join(clusteringPath, filename)
            ## Concatenate all new columns to the original DataFrame
            df_frac_combined = pd.concat([sampled_df] + clustering_columns, axis=1)
            df_frac_combined.to_csv(clusteringfile)
    
            dict_df_vScores[distType] = vScore_df
            dict_df_homogeneityScores[distType] = homogeneityScore_df
            dict_df_completenessScores[distType] = completenessScore_df
        
        ## Collect all scores for each cancer type
        dict_allScores_hierarchical[cancer_type] = [
                                                    dict_df_vScores,
                                                    dict_df_homogeneityScores,
                                                    dict_df_completenessScores,
                                                   ]
        
        ### Storing the files
        with open(f'TCGA_{cancer_type}_Clustering_Scores_Hierarchical_complete_{distType}_round{round}_exper01.pickle', 'wb') as handle:
            pickle.dump(dict_allScores_hierarchical, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    ## The end time
    end_time = time.time()
    print("Elapsed time:", (end_time - start_time)/60, f"minutes for {cancer_type}")
