"""
swap_gene_names.py 

2020-07-31
Updated 2020-10-29: moved mapping to a CSV in the db folder.

Adaptive MIRA Release 2 to tcrdist3 alphabeta_gammadelta_db.tsv
Generated 2020-07-31 Based on universal.py 
"""
import pandas as pd
from os.path import join as opj

from .paths import path_to_db

__all__ = ['adaptive_to_imgt']

mapping = pd.read_csv(opj(path_to_db, 'adaptive_imgt_mapping.csv'))

hum = mapping.loc[mapping['species'] == 'human'].set_index('adaptive')['imgt'].fillna('NA').to_dict()
mus = mapping.loc[mapping['species'] == 'mouse'].set_index('adaptive')['imgt'].fillna('NA').to_dict()

adaptive_to_imgt = {'human':hum,
                    'mouse':mus}

# This quickly makes a dictionary from  'TRGV7*00' -> 'TRGV7*01'
vdj00_to_imgt = dict()
vdj00_to_imgt['human']=\
	{f"{v.split('*')[0]}":v for k,v in adaptive_to_imgt['human'].items() }
vdj00_to_imgt['mouse']=\
	{f"{v.split('*')[0]}":v for k,v in adaptive_to_imgt['mouse'].items() }


