# emd_distance.py

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from pyemd import emd
import time

class EMDDistanceCalculator:
    def __init__(self, exposure_data, distance_matrix):
        """
        Initialize the EMDDistanceCalculator with exposure data and a signature distance matrix.
        
        Parameters:
        - exposure_data: Pandas DataFrame where rows are samples and columns are signature exposures.
        - distance_matrix: Pandas DataFrame containing the distance matrix between signatures.
        """
        self.exposure_data = exposure_data
        self.distance_matrix = distance_matrix

    def compute_emd(self, metric='cosine'):
        """
        Compute the EMD distance matrix for the given metric.
        
        Parameters:
        - metric: The type of distance to use ('cosine', 'hybrid', 'etiological').
        
        Returns:
        - dist_df: Pandas DataFrame containing the EMD distance matrix between samples.
        """
        if metric == 'cosine':
            SigDistMat = self.distance_matrix['cosine']
        elif metric == 'hybrid':
            SigDistMat = self.distance_matrix['hybrid']
        elif metric == 'etiological':
            SigDistMat = self.distance_matrix['etiological']
        else:
            raise ValueError("Metric must be 'cosine', 'hybrid', or 'etiological'")
        
        def emd_metric(p, q):
            return emd(np.ascontiguousarray(p, dtype=np.float64),
                       np.ascontiguousarray(q, dtype=np.float64),
                       np.ascontiguousarray(SigDistMat, dtype=np.float64))

        # Calculate pairwise EMD distances between samples
        dists = pdist(self.exposure_data, metric=emd_metric)
        dist_mat = squareform(dists)
        dist_df = pd.DataFrame(dist_mat, index=self.exposure_data.index, columns=self.exposure_data.index)

        return dist_df

def calculate_all_distances(exposure_data_dict, distance_matrix_dict, cancer_types):
    """
    Calculate EMD distance matrices for multiple cancer types and metrics.
    
    Parameters:
    - exposure_data_dict: Dictionary of DataFrames with cancer types as keys and exposure data as values.
    - distance_matrix_dict: Dictionary of distance matrices with metrics as keys and DataFrames as values.
    - cancer_types: List of cancer types to process.
    
    Returns:
    - dict_distDFs: Dictionary with cancer types as keys and lists of EMD distance DataFrames as values.
    """
    dict_distDFs = {}
    start_time = time.time()

    for cancer_type in cancer_types:
        distDFs = []
        exposure_data = exposure_data_dict[cancer_type]

        # Initialize the EMD calculator for each type of distance
        emd_calculator_cosine = EMDDistanceCalculator(exposure_data, {'cosine': distance_matrix_dict['cosine']})
        emd_calculator_hybrid = EMDDistanceCalculator(exposure_data, {'hybrid': distance_matrix_dict['hybrid']})
        emd_calculator_etiological = EMDDistanceCalculator(exposure_data, {'etiological': distance_matrix_dict['etiological']})

        # Calculate distances using different metrics
        distDFs.append(emd_calculator_cosine.compute_emd(metric='cosine'))
        distDFs.append(emd_calculator_hybrid.compute_emd(metric='hybrid'))
        distDFs.append(emd_calculator_etiological.compute_emd(metric='etiological'))

        dict_distDFs[cancer_type] = distDFs

        end_time = time.time()
        print("Elapsed time:", (end_time - start_time)/60, f"minutes for {cancer_type}")

    return dict_distDFs