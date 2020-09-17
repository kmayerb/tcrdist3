```python
import itertools
import math
import numpy as np
import os
import pandas as pd
import sys

from collections import Counter
from tcrdist.pgen import OlgaModel
from tcrsampler.sampler import TCRsampler
from tcrdist.repertoire import TCRrep
from tcrdist.automate import auto_pgen
from tcrdist.neighbors import compute_ecdf, bkgd_cntl_nn2
from tcrdist.background import make_gene_usage_counter, make_vj_matched_background, make_flat_vj_background, get_gene_frequencies, calculate_adjustment

ts = TCRsampler(default_background = 'britanova_human_beta_t_cb.tsv.sampler.tsv')
ts.build_background(stratify_by_subject = True, use_frequency = False)


# Functions for making a background set
ix =[['TRBV19*01', 'TRBJ2-5*01', 3],
['TRBV24-1*01', 'TRBJ2-4*01', 3],
['TRBV25-1*01', 'TRBJ2-4*01', 3],
['TRBV30*01', 'TRBJ2-3*01', 2],
['TRBV5-4*01', 'TRBJ2-3*01', 2],
['TRBV11-2*01', 'TRBJ2-2*01', 2],
['TRBV2*01', 'TRBJ1-5*01', 1],
['TRBV12-5*01', 'TRBJ2-7*01', 1],
['TRBV4-1*01', 'TRBJ1-6*01', 1],
['TRBV6-5*01', 'TRBJ1-6*01', 1],
['TRBV13*01', 'TRBJ2-3*01', 1],
['TRBV18*01', 'TRBJ2-3*01', 1],
['TRBV14*01', 'TRBJ2-7*01', 1],
['TRBV6-6*01', 'TRBJ2-7*01', 1],
['TRBV10-3*01', 'TRBJ2-3*01', 1],
['TRBV7-2*01', 'TRBJ2-1*01', 1],
['TRBV5-1*01', 'TRBJ2-1*01', 1]]

# For testing we  define a set of relatively rare TCRs 
flatten = lambda l: [item for sublist in l for item in sublist]
df_rare= pd.concat([pd.DataFrame({'cdr3_b_aa' : flatten(ts.sample([[x[0], x[1], x[2]]])) , 'v_b_gene':x[0], 'j_b_gene':x[1]}) for x in ix]).reset_index(drop = True)

# Make a VJ-Matched Background 
gene_usage_counter = make_gene_usage_counter(df_rare)
df_vj_bkgd = make_vj_matched_background(ts = ts,
	gene_usage_counter = gene_usage_counter, 	
	size = 100000, 
	recomb_type="VDJ", 
	chain_folder = "human_T_beta",
	cols = ['v_b_gene', 'j_b_gene', 'cdr3_b_aa'])
# Make a VJ flat 100K background
df_vj_flat = make_flat_vj_background(ts = ts, size = 100000)
df_bkgd = pd.concat([df_vj_bkgd, df_vj_flat]).reset_index(drop = True)
df_bkgd['weights'] = calculate_adjustment(df = df_bkgd)
```

```
         v_b_gene    j_b_gene             cdr3_b_aa        pV        pJ       pVJ   weights
0       TRBV19*01  TRBJ2-5*01      CASSIDLPGLEETQYF  0.006788  0.089505  0.000394  0.006760
1       TRBV19*01  TRBJ2-5*01      CASSIVGALILETQYF  0.006788  0.089505  0.000394  0.006760
2       TRBV19*01  TRBJ2-5*01   CASSIRRTGHGRREETQYF  0.006788  0.089505  0.000394  0.006760
3       TRBV19*01  TRBJ2-5*01      CASSIGGTNWTETQYF  0.006788  0.089505  0.000394  0.006760
4       TRBV19*01  TRBJ2-5*01        CAPLLNLAAITQYF  0.006788  0.089505  0.000394  0.006760
...           ...         ...                   ...       ...       ...       ...       ...
199837   TRBV9*01  TRBJ2-7*01      CASRTKIGGDSYEQYF  0.018002  0.217610  0.003250  4.811239
199838   TRBV9*01  TRBJ2-7*01         CASSDNRSYEQYF  0.018002  0.217610  0.003250  4.811239
199839   TRBV9*01  TRBJ2-7*01      CASRLGRGAMGYEQYF  0.018002  0.217610  0.003250  4.811239
199840   TRBV9*01  TRBJ2-7*01     CASSVAPQVSGSYEQYF  0.018002  0.217610  0.003250  4.811239
199841   TRBV9*01  TRBJ2-7*01  CASSVSGKEAGLGGIYEQYF  0.018002  0.217610  0.003250  4.811239
```



```python
tr_rare = TCRrep(cell_df = df_rare, organism = 'human', chains = ['beta'], store_all_cdr= False, compute_distances = False)
auto_pgen(tr_rare)
# Let's see the proportion of TCRs in a real massive dataset < 30 TCRDIST
df_check = ts.ref_df.rename(columns = {'v_reps':'v_b_gene','j_reps' :'j_b_gene', 'cdr3' : 'cdr3_b_aa'})[['v_b_gene','j_b_gene', 'cdr3_b_aa']].copy()
tr_check  = TCRrep(cell_df = df_check, organism = 'human', chains =['beta'], compute_distances = False, store_all_cdr = False)
tr_rare = TCRrep(cell_df = df_rare, organism = 'human', chains = ['beta'], store_all_cdr= False)
tr_rare.compute_rect_distances(df = tr_rare.clone_df, df2 = tr_check.clone_df )

tr_rare.clone_df['actual_hit_freq'] = np.sum(tr_rare.rw_beta <=30, axis = 1) / tr_rare.rw_beta.shape[1]
tr_rare.clone_df['actual_hit_count']= np.sum(tr_rare.rw_beta <=30, axis = 1)
df_save = tr_rare.clone_df.copy()

tr_bkgd = TCRrep(cell_df = df_bkgd, organism = 'human', chains =['beta'], compute_distances = False, store_all_cdr = False)
tr_rare.compute_rect_distances(df = tr_rare.clone_df, df2 = tr_bkgd.clone_df)
tr_rare.clone_df['raw_hit_freq'] = np.sum(tr_rare.rw_beta <=30, axis = 1) / tr_rare.rw_beta.shape[1]
tr_rare.clone_df['raw_hit_count']= np.sum(tr_rare.rw_beta <=30, axis = 1)
tr_rare.clone_df['raw_bias'] = tr_rare.clone_df['raw_hit_freq'] / tr_rare.clone_df['actual_hit_freq'] 
weights = tr_bkgd.clone_df.weights
tr_rare.clone_df['adj_hit_count'] =  np.sum((1*(tr_rare.rw_beta <=30) * weights[None,:]), axis = 1)
tr_rare.clone_df['adj_hit_freq'] =  np.sum((1*(tr_rare.rw_beta <=30) * weights[None,:]), axis = 1) / tr_rare.rw_beta.shape[1]
tr_rare.clone_df['adj_bias'] = tr_rare.clone_df['adj_hit_freq'] / tr_rare.clone_df['actual_hit_freq'] 


auto_pgen(tr_rare)
df_save2 = tr_rare.clone_df.copy()
df_save2['actual_hit_freq'] = df_save2['actual_hit_freq'].apply(np.format_float_scientific,precision =1)
df_save2['raw_hit_freq'] = df_save2['raw_hit_freq'].apply(np.format_float_scientific,precision =1)
df_save2['adj_hit_freq'] = df_save2['adj_hit_freq'].apply(np.format_float_scientific,precision =1)
df_save2[['cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'pgen_cdr3_b_aa', 'raw_hit_count', 'actual_hit_freq','raw_hit_freq','adj_hit_freq','raw_bias', 'adj_bias' ]]

```

## flat-vj + vj-matched

### Using a Composite Background 100K VJ Matched Set + 100K Flat Set Works Well


```
             cdr3_b_aa     v_b_gene    j_b_gene  pgen_cdr3_b_aa  raw_hit_count actual_hit_freq raw_hit_freq adj_hit_freq   raw_bias  adj_bias
0       CAISESANTDTQYF  TRBV10-3*01  TRBJ2-3*01    6.915652e-10             24         1.2e-04      1.2e-04      7.1e-05   0.968908  0.569920
1       CARTSGRETDTQYF   TRBV5-4*01  TRBJ2-3*01    1.636423e-09             21         4.5e-06      1.1e-04      1.3e-06  23.644040  0.294601
2       CASGPGQGPYEQYF  TRBV12-5*01  TRBJ2-7*01    2.254880e-08             25         2.3e-05      1.3e-04      5.1e-06   5.507152  0.222898
3         CASHRDGQPQHF     TRBV2*01  TRBJ1-5*01    4.744571e-09             17         9.9e-06      8.6e-05      2.7e-06   8.613186  0.272006
4        CASRQGSYNEQFF   TRBV5-1*01  TRBJ2-1*01    2.865843e-07             37         1.2e-04      1.9e-04      2.3e-04   1.524093  1.852166
5      CASSAGQVNTGELFF  TRBV11-2*01  TRBJ2-2*01    1.127123e-08             22         3.7e-06      1.1e-04      1.5e-06  29.723935  0.407522
6        CASSDRAKNIQYF  TRBV25-1*01  TRBJ2-4*01    2.769610e-07            333         3.9e-05      1.7e-03      1.3e-05  42.713193  0.320593
7       CASSESLSENIQYF  TRBV25-1*01  TRBJ2-4*01    3.299466e-10             45         6.5e-06      2.3e-04      1.7e-06  35.076322  0.263273
8   CASSHGLGGFPFTDTQYF   TRBV5-4*01  TRBJ2-3*01    2.774562e-13              0         2.5e-07       0.e+00       0.e+00   0.000000  0.000000
9          CASSIQETQYF    TRBV19*01  TRBJ2-5*01    3.226682e-07            550         1.1e-04      2.8e-03      4.4e-05  24.390538  0.389876
10     CASSISHRGLETQYF    TRBV19*01  TRBJ2-5*01    3.145629e-10             20         3.7e-06      1.0e-04      6.8e-07  27.021759  0.182678
11       CASSLEYTGELFF  TRBV11-2*01  TRBJ2-2*01    6.950468e-08            179         3.9e-05      9.0e-04      1.3e-05  22.815542  0.327919
12      CASSLIDGIGEQYF   TRBV6-6*01  TRBJ2-7*01    3.803834e-11              0          4.e-06       0.e+00       0.e+00   0.000000  0.000000
13      CASSPGTGGSIQYF  TRBV25-1*01  TRBJ2-4*01    1.346750e-09              2         1.7e-06      1.0e-05      7.6e-08   5.790377  0.043461
14    CASSPLTGAGTDTQYF    TRBV18*01  TRBJ2-3*01    2.255192e-09             23         3.2e-05      1.2e-04      1.1e-05   3.670278  0.352780
15    CASSPRRRGGADTQYF    TRBV13*01  TRBJ2-3*01    6.799742e-10              1         9.9e-07      5.0e-06      3.8e-07   5.066580  0.379799
16       CASSPVGLYEQYF    TRBV14*01  TRBJ2-7*01    5.567815e-08             32         5.2e-05      1.6e-04      5.0e-05   3.088201  0.959476
17        CASSQTPETQYF    TRBV19*01  TRBJ2-5*01    5.264363e-08            104         2.2e-05      5.2e-04      9.9e-06  24.226405  0.459236
18     CASSRDGRRNSPLHF   TRBV4-1*01  TRBJ1-6*01    6.167066e-10              0         8.7e-06       0.e+00       0.e+00   0.000000  0.000000
19     CASSWTSGRANEQFF   TRBV7-2*01  TRBJ2-1*01    1.009388e-08              1         6.3e-05      5.0e-06      3.1e-06   0.079476  0.049637
20      CASSYRGPDSPLHF   TRBV6-5*01  TRBJ1-6*01    9.649681e-09              0         4.7e-05       0.e+00       0.e+00   0.000000  0.000000
21      CATGEKIAKNIQYF  TRBV24-1*01  TRBJ2-4*01    1.041126e-10             71         4.2e-06      3.6e-04      2.6e-06  84.641688  0.609329
22      CATSALLATNIQYF  TRBV24-1*01  TRBJ2-4*01    4.811646e-12             14         1.5e-06      7.0e-05      5.1e-07  47.288079  0.340423
23     CATSDRDRAKNIQYF  TRBV24-1*01  TRBJ2-4*01    6.159715e-09            151         9.9e-06      7.6e-04      5.5e-06  76.505356  0.550756
24       CAWSPTLADTQYF    TRBV30*01  TRBJ2-3*01    4.336671e-10             19         1.2e-05      9.6e-05      1.1e-06   7.858369  0.094562
25   CAWTPLAGGRFTDTQYF    TRBV30*01  TRBJ2-3*01    9.201160e-12              1         2.5e-07      5.0e-06      6.1e-08  20.266320  0.243869
```

## flat-vj + vj-matched + real

### Adding and additional 100K Real Seqs has virtually no effect on bias


```python

df_real = ts.ref_df.sample(100000, random_state = 1).rename(columns = {'v_reps':'v_b_gene','j_reps' :'j_b_gene', 'cdr3' : 'cdr3_b_aa'})[['v_b_gene','j_b_gene', 'cdr3_b_aa']]
df_real = get_gene_frequencies(ts = ts, df = df_real) 
df_real['weights'] =  calculate_adjustment(df = df_real)
df_bkgd2 = pd.concat([df_bkgd, df_real]).reset_index(drop = True)

tr_bkgd = TCRrep(cell_df = df_bkgd2, organism = 'human', chains =['beta'], compute_distances = False, store_all_cdr = False)
tr_rare.compute_rect_distances(df = tr_rare.clone_df, df2 = tr_bkgd.clone_df)
tr_rare.clone_df['raw_hit_freq'] = np.sum(tr_rare.rw_beta <=30, axis = 1) / tr_rare.rw_beta.shape[1]
tr_rare.clone_df['raw_hit_count']= np.sum(tr_rare.rw_beta <=30, axis = 1)
tr_rare.clone_df['raw_bias'] = tr_rare.clone_df['raw_hit_freq'] / tr_rare.clone_df['actual_hit_freq'] 
weights = tr_bkgd.clone_df.weights
tr_rare.clone_df['adj_hit_count'] =  np.sum((1*(tr_rare.rw_beta <=30) * weights[None,:]), axis = 1)
tr_rare.clone_df['adj_hit_freq'] =  np.sum((1*(tr_rare.rw_beta <=30) * weights[None,:]), axis = 1) / tr_rare.rw_beta.shape[1]
tr_rare.clone_df['adj_bias'] = tr_rare.clone_df['adj_hit_freq'] / tr_rare.clone_df['actual_hit_freq'] 
auto_pgen(tr_rare)
df_save3 = tr_rare.clone_df.copy()
df_save3['actual_hit_freq'] = df_save3['actual_hit_freq'].apply(np.format_float_scientific,precision =1)
df_save3['raw_hit_freq'] = df_save3['raw_hit_freq'].apply(np.format_float_scientific,precision =1)
df_save3['adj_hit_freq'] = df_save3['adj_hit_freq'].apply(np.format_float_scientific,precision =1)
df_save3[['cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'pgen_cdr3_b_aa', 'raw_hit_count', 'actual_hit_freq','raw_hit_freq','adj_hit_freq','raw_bias', 'adj_bias' ]]
```

```
             cdr3_b_aa     v_b_gene    j_b_gene  pgen_cdr3_b_aa  raw_hit_count actual_hit_freq raw_hit_freq adj_hit_freq   raw_bias  adj_bias
0       CAISESANTDTQYF  TRBV10-3*01  TRBJ2-3*01    6.915652e-10             25         1.2e-04      1.2e-04      7.5e-05   0.998466  0.600501
1       CARTSGRETDTQYF   TRBV5-4*01  TRBJ2-3*01    1.636423e-09             21         4.5e-06      1.0e-04      1.3e-06  23.390731  0.291445
2       CASGPGQGPYEQYF  TRBV12-5*01  TRBJ2-7*01    2.254880e-08             25         2.3e-05      1.2e-04      5.0e-06   5.448152  0.220510
3         CASHRDGQPQHF     TRBV2*01  TRBJ1-5*01    4.744571e-09             17         9.9e-06      8.5e-05      2.7e-06   8.520909  0.269092
4        CASRQGSYNEQFF   TRBV5-1*01  TRBJ2-1*01    2.865843e-07             37         1.2e-04      1.8e-04      2.2e-04   1.507765  1.832323
5      CASSAGQVNTGELFF  TRBV11-2*01  TRBJ2-2*01    1.127123e-08             22         3.7e-06      1.1e-04      1.5e-06  29.405491  0.403156
6        CASSDRAKNIQYF  TRBV25-1*01  TRBJ2-4*01    2.769610e-07            333         3.9e-05      1.7e-03      1.2e-05  42.255588  0.317158
7       CASSESLSENIQYF  TRBV25-1*01  TRBJ2-4*01    3.299466e-10             45         6.5e-06      2.2e-04      1.7e-06  34.700535  0.260452
8   CASSHGLGGFPFTDTQYF   TRBV5-4*01  TRBJ2-3*01    2.774562e-13              0         2.5e-07       0.e+00       0.e+00   0.000000  0.000000
9          CASSIQETQYF    TRBV19*01  TRBJ2-5*01    3.226682e-07            551         1.1e-04      2.7e-03      4.8e-05  24.173103  0.426102
10     CASSISHRGLETQYF    TRBV19*01  TRBJ2-5*01    3.145629e-10             20         3.7e-06       1.e-04      6.7e-07  26.732264  0.180721
11       CASSLEYTGELFF  TRBV11-2*01  TRBJ2-2*01    6.950468e-08            179         3.9e-05      8.9e-04      1.3e-05  22.571110  0.324406
12      CASSLIDGIGEQYF   TRBV6-6*01  TRBJ2-7*01    3.803834e-11              0          4.e-06       0.e+00       0.e+00   0.000000  0.000000
13      CASSPGTGGSIQYF  TRBV25-1*01  TRBJ2-4*01    1.346750e-09              2         1.7e-06       1.e-05      7.5e-08   5.728342  0.042995
14    CASSPLTGAGTDTQYF    TRBV18*01  TRBJ2-3*01    2.255192e-09             23         3.2e-05      1.1e-04      1.1e-05   3.630957  0.349000
15    CASSPRRRGGADTQYF    TRBV13*01  TRBJ2-3*01    6.799742e-10              1         9.9e-07       5.e-06      3.7e-07   5.012300  0.375730
16       CASSPVGLYEQYF    TRBV14*01  TRBJ2-7*01    5.567815e-08             32         5.2e-05      1.6e-04      4.9e-05   3.055116  0.949197
17        CASSQTPETQYF    TRBV19*01  TRBJ2-5*01    5.264363e-08            104         2.2e-05      5.2e-04      9.8e-06  23.966858  0.454316
18     CASSRDGRRNSPLHF   TRBV4-1*01  TRBJ1-6*01    6.167066e-10              0         8.7e-06       0.e+00       0.e+00   0.000000  0.000000
19     CASSWTSGRANEQFF   TRBV7-2*01  TRBJ2-1*01    1.009388e-08              1         6.3e-05       5.e-06      3.1e-06   0.078624  0.049105
20      CASSYRGPDSPLHF   TRBV6-5*01  TRBJ1-6*01    9.649681e-09              0         4.7e-05       0.e+00       0.e+00   0.000000  0.000000
21      CATGEKIAKNIQYF  TRBV24-1*01  TRBJ2-4*01    1.041126e-10             71         4.2e-06      3.5e-04      2.5e-06  83.734886  0.602801
22      CATSALLATNIQYF  TRBV24-1*01  TRBJ2-4*01    4.811646e-12             14         1.5e-06       7.e-05      5.0e-07  46.781462  0.336776
23     CATSDRDRAKNIQYF  TRBV24-1*01  TRBJ2-4*01    6.159715e-09            151         9.9e-06      7.5e-04      5.4e-06  75.685723  0.544856
24       CAWSPTLADTQYF    TRBV30*01  TRBJ2-3*01    4.336671e-10             19         1.2e-05      9.5e-05      1.1e-06   7.774179  0.093548
25   CAWTPLAGGRFTDTQYF    TRBV30*01  TRBJ2-3*01    9.201160e-12              1         2.5e-07       5.e-06       6.e-08  20.049198  0.241257
```

## flat-vj only 

### But, using only the 100K flat background has too many zeros.

```
df_bkgd0 = pd.concat([df_vj_flat]).reset_index(drop = True)
df_bkgd0['weights'] = calculate_adjustment(df = df_bkgd0)

tr_bkgd = TCRrep(cell_df = df_bkgd0, organism = 'human', chains =['beta'], compute_distances = False, store_all_cdr = False)
tr_rare.compute_rect_distances(df = tr_rare.clone_df, df2 = tr_bkgd.clone_df)
tr_rare.clone_df['raw_hit_freq'] = np.sum(tr_rare.rw_beta <=30, axis = 1) / tr_rare.rw_beta.shape[1]
tr_rare.clone_df['raw_hit_count']= np.sum(tr_rare.rw_beta <=30, axis = 1)
tr_rare.clone_df['raw_bias'] = tr_rare.clone_df['raw_hit_freq'] / tr_rare.clone_df['actual_hit_freq'] 
weights = tr_bkgd.clone_df.weights
tr_rare.clone_df['adj_hit_count'] =  np.sum((1*(tr_rare.rw_beta <=30) * weights[None,:]), axis = 1)
tr_rare.clone_df['adj_hit_freq'] =  np.sum((1*(tr_rare.rw_beta <=30) * weights[None,:]), axis = 1) / tr_rare.rw_beta.shape[1]
tr_rare.clone_df['adj_bias'] = tr_rare.clone_df['adj_hit_freq'] / tr_rare.clone_df['actual_hit_freq'] 


auto_pgen(tr_rare)
df_save3 = tr_rare.clone_df.copy()
df_save3['actual_hit_freq'] = df_save3['actual_hit_freq'].apply(np.format_float_scientific,precision =1)
df_save3['raw_hit_freq'] = df_save3['raw_hit_freq'].apply(np.format_float_scientific,precision =1)
df_save3['adj_hit_freq'] = df_save3['adj_hit_freq'].apply(np.format_float_scientific,precision =1)
df_save3[['cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'pgen_cdr3_b_aa', 'raw_hit_count', 'actual_hit_freq','raw_hit_freq','adj_hit_freq','raw_bias', 'adj_bias' ]]
```

```
             cdr3_b_aa     v_b_gene    j_b_gene  pgen_cdr3_b_aa  raw_hit_count actual_hit_freq raw_hit_freq adj_hit_freq  raw_bias  adj_bias
0       CAISESANTDTQYF  TRBV10-3*01  TRBJ2-3*01    6.915652e-10              2         1.2e-04       2.e-05      8.2e-05  0.160519  0.660710
1       CARTSGRETDTQYF   TRBV5-4*01  TRBJ2-3*01    1.636423e-09              0         4.5e-06       0.e+00       0.e+00  0.000000  0.000000
2       CASGPGQGPYEQYF  TRBV12-5*01  TRBJ2-7*01    2.254880e-08              2         2.3e-05       2.e-05      1.2e-05  0.875875  0.522636
3         CASHRDGQPQHF     TRBV2*01  TRBJ1-5*01    4.744571e-09              0         9.9e-06       0.e+00       0.e+00  0.000000  0.000000
4        CASRQGSYNEQFF   TRBV5-1*01  TRBJ2-1*01    2.865843e-07              2         1.2e-04       2.e-05      3.6e-04  0.163782  2.932869
5      CASSAGQVNTGELFF  TRBV11-2*01  TRBJ2-2*01    1.127123e-08              1         3.7e-06       1.e-05       4.e-06  2.686017  1.067122
6        CASSDRAKNIQYF  TRBV25-1*01  TRBJ2-4*01    2.769610e-07              4         3.9e-05       4.e-05      1.3e-05  1.020007  0.330629
7       CASSESLSENIQYF  TRBV25-1*01  TRBJ2-4*01    3.299466e-10              1         6.5e-06       1.e-05      3.2e-06  1.549625  0.502301
8   CASSHGLGGFPFTDTQYF   TRBV5-4*01  TRBJ2-3*01    2.774562e-13              0         2.5e-07       0.e+00       0.e+00  0.000000  0.000000
9          CASSIQETQYF    TRBV19*01  TRBJ2-5*01    3.226682e-07              7         1.1e-04       7.e-05      3.4e-05  0.617137  0.302268
10     CASSISHRGLETQYF    TRBV19*01  TRBJ2-5*01    3.145629e-10              0         3.7e-06       0.e+00       0.e+00  0.000000  0.000000
11       CASSLEYTGELFF  TRBV11-2*01  TRBJ2-2*01    6.950468e-08              3         3.9e-05       3.e-05      8.6e-06  0.760194  0.218122
12      CASSLIDGIGEQYF   TRBV6-6*01  TRBJ2-7*01    3.803834e-11              0          4.e-06       0.e+00       0.e+00  0.000000  0.000000
13      CASSPGTGGSIQYF  TRBV25-1*01  TRBJ2-4*01    1.346750e-09              0         1.7e-06       0.e+00       0.e+00  0.000000  0.000000
14    CASSPLTGAGTDTQYF    TRBV18*01  TRBJ2-3*01    2.255192e-09              0         3.2e-05       0.e+00       0.e+00  0.000000  0.000000
15    CASSPRRRGGADTQYF    TRBV13*01  TRBJ2-3*01    6.799742e-10              0         9.9e-07       0.e+00       0.e+00  0.000000  0.000000
16       CASSPVGLYEQYF    TRBV14*01  TRBJ2-7*01    5.567815e-08              1         5.2e-05       1.e-05      3.0e-05  0.191858  0.580979
17        CASSQTPETQYF    TRBV19*01  TRBJ2-5*01    5.264363e-08              3         2.2e-05       3.e-05      1.2e-05  1.389319  0.565865
18     CASSRDGRRNSPLHF   TRBV4-1*01  TRBJ1-6*01    6.167066e-10              0         8.7e-06       0.e+00       0.e+00  0.000000  0.000000
19     CASSWTSGRANEQFF   TRBV7-2*01  TRBJ2-1*01    1.009388e-08              0         6.3e-05       0.e+00       0.e+00  0.000000  0.000000
20      CASSYRGPDSPLHF   TRBV6-5*01  TRBJ1-6*01    9.649681e-09              0         4.7e-05       0.e+00       0.e+00  0.000000  0.000000
21      CATGEKIAKNIQYF  TRBV24-1*01  TRBJ2-4*01    1.041126e-10              0         4.2e-06       0.e+00       0.e+00  0.000000  0.000000
22      CATSALLATNIQYF  TRBV24-1*01  TRBJ2-4*01    4.811646e-12              0         1.5e-06       0.e+00       0.e+00  0.000000  0.000000
23     CATSDRDRAKNIQYF  TRBV24-1*01  TRBJ2-4*01    6.159715e-09              2         9.9e-06       2.e-05      6.2e-06  2.014513  0.626624
24       CAWSPTLADTQYF    TRBV30*01  TRBJ2-3*01    4.336671e-10              1         1.2e-05       1.e-05      3.5e-06  0.822250  0.286786
25   CAWTPLAGGRFTDTQYF    TRBV30*01  TRBJ2-3*01    9.201160e-12              0         2.5e-07       0.e+00       0.e+00  0.000000  0.000000
```
