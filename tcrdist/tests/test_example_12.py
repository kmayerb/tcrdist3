import pandas as pd 
from tcrdist.repertoire import TCRrep

import sys
import numpy as np
import pytest

pytestmark = pytest.mark.skipif(
    sys.version_info >= (3, 11) and int(np.__version__.split(".", 1)[0]) >= 2,
    reason="Fails on Python 3.11 with numpy >= 2.0 (uses private API)"
)

def test_example_12():
	"""
	Illustrate how to save a reload from and archive file
	"""
	import pandas as pd
	from tcrdist.repertoire import TCRrep
	import numpy as np

	df = pd.read_csv("dash.csv").head(100)
	tr = TCRrep(cell_df = df,
		organism = 'mouse',
		chains = ['alpha','beta'],
		db_file = 'alphabeta_gammadelta_db.tsv',
		store_all_cdr=False,
		archive_result=True,
		archive_name = "example_archive")

	tr2 = TCRrep(cell_df = None,
		organism = 'mouse',
		chains = ['alpha','beta'],
		db_file = 'alphabeta_gammadelta_db.tsv',
		blank = True,
		archive_name = "example_archive")
	tr2.rebuild()

	# Check that all atrributes are the same after rebuild, except metrics which can't be zipped
	for k in tr2.__dict__.keys():
		print(k)
		if k in ['all_genes','metrics_a','metrics_b','metrics_d', 'metrics_g',
				'kargs_a','kargs_b','kargs_d','kargs_g']:
			pass
		else:
			assert np.all(getattr(tr, k) == getattr(tr2, k) )

	for k in ['all_genes','metrics_a','metrics_b', 'kargs_a','kargs_b']:
		assert isinstance(getattr(tr2, k), dict)
		assert set(getattr(tr, k).keys()) - set(getattr(tr2, k).keys()) == set()
