import pandas as pd
from tcrdist.repertoire import TCRrep

import pandas as pd
from tcrdist.repertoire import TCRrep

df = pd.read_csv("dash.csv")
tr = TCRrep(cell_df = df,
            organism = 'mouse',
            chains = ['alpha','beta'],
            db_file = 'alphabeta_gammadelta_db.tsv',
            compute_distances = False)

tr.cpus = 2
tr.compute_sparse_rect_distances(radius = 50, chunk_size = 100)
tr.rw_beta
"""<1920x1920 sparse matrix of type '<class 'numpy.int16'>'
	with 108846 stored elements in Compressed Sparse Row format>
"""
print(tr.rw_beta)
"""
  (0, 0)  -1
  (1, 1)  -1
  (1, 470)  24
  (1, 472)  24
  (2, 2)  -1
  : :
  (1919, 1911)  24
  (1919, 1912)  38
  (1919, 1918)  12
  (1919, 1919)  -1
"""