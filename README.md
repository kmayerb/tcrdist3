# tcrdist3

![Python application](https://github.com/kmayerb/tcrdist3/workflows/Python%20application/badge.svg?event=push) [![Coverage Status](https://coveralls.io/repos/github/kmayerb/tcrdist3/badge.svg?branch=master)](https://coveralls.io/github/kmayerb/tcrdist3?branch=master)[![Documentation Status](https://readthedocs.org/projects/tcrdist3/badge/?version=latest)](https://tcrdist3.readthedocs.io/en/latest/?badge=latest)
[![Docker Repository on Quay](https://quay.io/repository/kmayerb/tcrdist3/status "Docker Repository on Quay")](https://quay.io/repository/kmayerb/tcrdist3)

Flexible distance measures for comparing T cell receptors 

tcrdist3 is a python API-enabled toolkit for analyzing T-cell receptor repertoires. Some of the functionality and code is adapted from the original tcr-dist package which was released with the publication of Dash et al. Nature (2017) doi:10.1038/nature22383. This package contains a new API for computing tcrdistance measures as well as new features for biomarker development ([bioRxiv (2020)](https://www.biorxiv.org/content/10.1101/2020.12.24.424260v1)). The package has been expanded to include gamma-delta TCRs; it has also been recoded to increase CPU efficiency using numba, a high-performance just-in-time compiler.

<img src="https://user-images.githubusercontent.com/46639063/103338268-aa3ee180-4a32-11eb-8149-056fb385b33b.gif" width="720">

## Installation

[![PyPI version](https://badge.fury.io/py/tcrdist3.svg)](https://badge.fury.io/py/tcrdist3)

```
pip install tcrdist3
```

or 

```
pip install git+https://github.com/kmayerb/tcrdist3.git@0.2.2
```

## Docker
[![Docker Repository on Quay](https://quay.io/repository/kmayerb/tcrdist3/status "Docker Repository on Quay")](https://quay.io/repository/kmayerb/tcrdist3)

```
docker pull quay.io/kmayerb/tcrdist3:0.2.2
```
## User-Contributed Colab Notebook Examples Using tcrdist3

### 1. Example K Nearest Neighbor Classification using tcrdist3 

[![open in colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1boqbGZjJqt_di3-3ygatHO-G4L7t4rY-?usp=sharing) (Author: Liel Cohen-Lavi). This notebook illustrates how to integrate tcrdist3 with scikit-learn's implementation of K Nearest Neighbor classification. TCRdist-based KNN classification performance on a set of labeled receptors is assessed with cross-validation or training/test splits   This simple method is proposed as a quickly implementable benchmark for the performance of more computationally intensive TCR-epitope specificity prediction approaches. 

## Package Documentation
[![Documentation Status](https://readthedocs.org/projects/tcrdist3/badge/?version=latest)](https://tcrdist3.readthedocs.io/en/latest/?badge=latest)

More documentation can be found at [tcrdist3.readthedocs](https://tcrdist3.readthedocs.io/).

## Basic Usage

```python
import pandas as pd
from tcrdist.repertoire import TCRrep

df = pd.read_csv("dash.csv")
tr = TCRrep(cell_df = df, 
            organism = 'mouse', 
            chains = ['alpha','beta'], 
            db_file = 'alphabeta_gammadelta_db.tsv')

tr.pw_alpha
tr.pw_beta
tr.pw_cdr3_a_aa
tr.pw_cdr3_b_aa

from tcrdist.public import _neighbors_fixed_radius
_neighbors_fixed_radius(tr.pw_beta, 50)         
```

### Sparse Matrix Representation 

```python
import pandas as pd
from tcrdist.repertoire import TCRrep
from tcrdist.breadth import get_safe_chunk

df = pd.read_csv("dash.csv")
tr = TCRrep(cell_df = df[['subject','epitope','count','v_b_gene','j_b_gene','cdr3_b_aa','cdr3_b_nucseq']], 
            organism = 'mouse', 
            chains = ['beta'], 
            compute_distances = False)

# Set to desired number of CPUs
tr.cpus = 2

# Identify a safe chunk size based on input data shape and target number of 
# pairwise distance to be temporarily held in memory per node. 
safe_chunk_size = get_safe_chunk(
            tr.clone_df.shape[0], 
            tr.clone_df.shape[0], 
            target = 10**7) 

tr.compute_sparse_rect_distances(
        df = tr.clone_df, 
        radius=50,
        chunk_size = safe_chunk_size)

print(tr.rw_beta)
```

## Citing

##### TCR meta-clonotypes for biomarker discovery with tcrdist3 enabled identification of public, HLA-restricted clusters of SARS-CoV-2 TCRs

Mayer-Blackwell K, Schattgen S, Cohen-Lavi L, Crawford JC, Souquette A, Gaevert JA, Hertz T, Thomas PG, Bradley PH, Fiore-Gartland A. [eLife (2021)](https://elifesciences.org/articles/68605).


##### Quantifiable predictive features define epitope-specific T cell receptor repertoires

Pradyot Dash, Andrew J. Fiore-Gartland, Tomer Hertz, George C. Wang, Shalini Sharma, Aisha Souquette, Jeremy Chase Crawford, E. Bridie Clemens, Thi H. O. Nguyen, Katherine Kedzierska, Nicole L. La Gruta, Philip Bradley & Paul G. Thomas [Nature (2017)](https://doi.org/10.1038/nature22383).
