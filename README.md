# 🧬 Mutational Signature Exposure Analysis with EMD

<p align="center">
  <img src="https://img.shields.io/badge/Python-3.11-blue?style=flat-square&logo=python"/>
  <img src="https://img.shields.io/badge/Status-Research-green?style=flat-square"/>
  <img src="https://img.shields.io/badge/License-MIT-yellow?style=flat-square"/>
</p>

## 📋 Overview

This project implements **Earth Mover's Distance (EMD)** based approaches for analyzing cancer mutational signature exposures. The key innovation is **hierarchical EMD (hEMD)** which incorporates biological etiology information to create more meaningful patient distance metrics.

### Key Features

- 🔬 **Etiology-aware distance metrics**: Incorporates mutational signature etiology (e.g., DNA repair deficiency, tobacco exposure) into distance calculations
- 📊 **Patient clustering**: Groups cancer patients based on their mutational signature profiles
- 🔄 **Multi-cohort analysis**: Supports TCGA and Hartwig datasets
- 📈 **Evaluation pipeline**: Comprehensive clustering evaluation with multiple metrics

## 🗂️ Repository Structure

```
Proj_MutSig_Exposure_EMD/
├── EMD_distMat_calculation.py          # Core EMD calculation module
├── example_usage_of_EMD_distMat_calculation.py  # Usage examples
├── _00_Prepare_*.ipynb                 # Data preparation notebooks
├── _01_Explore_data_*.ipynb            # Data exploration & distance calculation
├── _02_draw_clustermap_*.ipynb         # Visualization
├── _03_Demo_of_Clustering_*.py         # Clustering analysis
├── _04_evaluation_of_clustering_*.ipynb # Evaluation metrics
└── requirements.txt                    # Dependencies
```

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/rcwang2024/Proj_MutSig_Exposure_EMD.git
cd Proj_MutSig_Exposure_EMD

# Install dependencies
pip install -r requirements.txt
```

### Basic Usage

```python
import pandas as pd
from EMD_distMat_calculation import EMDDistanceCalculator, calculate_all_distances

# Patient x Signature exposure matrix (rows = samples, columns = SBS signatures)
exposures = pd.read_csv("breast_cancer_exposure.csv", index_col=0)

# Signature x Signature ground-distance matrices
sig_dist = {
    "cosine":      pd.read_csv("Cosine_distance_between_SBS_signatures.csv", index_col=0),
    "hybrid":      pd.read_csv("Hybrid_distance_between_SBS_signatures.csv", index_col=0),
    "etiological": pd.read_csv("Functional_distance_between_SBS_signatures.csv", index_col=0),
}

# Sample x Sample EMD distance matrix using the etiology-aware ground distance
emd_etio = EMDDistanceCalculator(exposures, sig_dist).compute_emd(metric="etiological")

# Or: all metrics for several cancer types at once
all_dists = calculate_all_distances({"Breast-cancer": exposures}, sig_dist, ["Breast-cancer"])
```

See [`example_usage_of_EMD_distMat_calculation.py`](example_usage_of_EMD_distMat_calculation.py) for a complete example.

## 📐 Methodology

### Earth Mover's Distance (EMD)

EMD measures the minimum "work" needed to transform one distribution into another, considering both the amount to move and the ground distance.

### Hierarchical EMD (hEMD)

Our approach uses a **hierarchical cost matrix** based on mutational signature etiology:

- **Same etiology**: Low cost (e.g., SBS4 ↔ SBS92, both tobacco-related)
- **Different etiology**: Higher cost reflecting biological difference

This enables biologically meaningful patient clustering beyond standard cosine distance.

## 📊 Analysis Pipeline

1. **Data Preparation** (`_00_*.ipynb`): Load signature exposures, prepare etiology mappings
2. **Distance Calculation** (`_01_*.ipynb`): Compute EMD distance matrices
3. **Visualization** (`_02_*.ipynb`): Generate cluster heatmaps
4. **Clustering** (`_03_*.py`): Hierarchical clustering analysis
5. **Evaluation** (`_04_*.ipynb`): Assess clustering quality with ARI, NMI, etc.

## 🛠️ Technologies

| Category | Tools |
|----------|-------|
| **Core** | Python 3.11, NumPy, SciPy |
| **Data** | Pandas |
| **ML** | Scikit-learn |
| **Visualization** | Matplotlib, Seaborn |
| **Optimal transport** | pyemd |

## 📚 References

- COSMIC Mutational Signatures: [https://cancer.sanger.ac.uk/signatures/](https://cancer.sanger.ac.uk/signatures/)
- pyemd: [https://github.com/wmayner/pyemd](https://github.com/wmayner/pyemd)
- Related project: [Distance Analysis of Mutational Signatures](https://github.com/rcwang2024/Proj_MutSig_distance_comparison) (ISMB/ECCB 2023 poster) — benchmark of signature-level distance metrics that motivated this work

## 📄 License

MIT License - see [LICENSE](LICENSE) for details.
