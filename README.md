# tcrdist3

![Python application](https://github.com/kmayerb/tcrdist3/workflows/Python%20application/badge.svg?event=push)

[![Coverage Status](https://coveralls.io/repos/github/kmayerb/tcrdist3/badge.svg?branch=master)](https://coveralls.io/github/kmayerb/tcrdist3?branch=master)


Flexible distance measures for comparing T cell receptors 


## Installation

```
pip install tcrdist3
```

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
```

## Citing

##### Quantifiable predictive features define epitope-specific T cell receptor repertoires

Pradyot Dash, Andrew J. Fiore-Gartland, Tomer Hertz, George C. Wang, Shalini Sharma, Aisha Souquette, Jeremy Chase Crawford, E. Bridie Clemens, Thi H. O. Nguyen, Katherine Kedzierska, Nicole L. La Gruta, Philip Bradley & Paul G. Thomas [Nature (2017)](https://doi.org/10.1038/nature22383).
