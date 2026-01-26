import pandas as pd
from EMD_distMat_calculation import calculate_all_distances

# Load the exposure data and distance matrices
dict_allFrac_ss = {
    'Breast-cancer': pd.read_csv('breast_cancer_exposure.csv'),
    'Lung-AdenoCa': pd.read_csv('lung_adenoca_exposure.csv'),
    # Add other cancer types as needed
}

distance_matrix_dict = {
    'cosine': pd.read_csv('Cosine_distance_between_SBS_signatures.csv', index_col=0),
    'hybrid': pd.read_csv('Hybrid_distance_between_SBS_signatures.csv', index_col=0),
    'etiological': pd.read_csv('Functional_distance_between_SBS_signatures.csv', index_col=0)
}

# Define the cancer types
cancer_types = ['Breast-cancer', 'Lung-AdenoCa']

# Calculate the EMD distance matrices
dict_distDFs = calculate_all_distances(dict_allFrac_ss, distance_matrix_dict, cancer_types)
