import pytest



def test_import_vdjtools_beta_w_validation():
    import pandas as pd
    import numpy as np
    import os
    from tcrdist.paths import path_to_base
    from tcrdist.vdjtools_funcs import import_vdjtools
    from tcrdist.repertoire import TCRrep

    # Reformat vdj_tools input format for tcrdist3
    vdj_tools_file_beta = os.path.join(path_to_base, 'tcrdist','data','formats','vdj.M_15_CD8_beta.clonotypes.TRB.txt.gz')
    df_beta = import_vdjtools(   vdj_tools_file = vdj_tools_file_beta ,
                        chain = 'beta',
                        organism = 'human',
                        db_file = 'alphabeta_gammadelta_db.tsv',
                        validate = True)
    assert np.all(df_beta.columns == ['count', 'freq', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'cdr3_b_nucseq','valid_v', 'valid_j', 'valid_cdr3'])
    
    # Can be directly imported into a TCRrep instance.
    tr = TCRrep(
        cell_df = df_beta[['count', 'freq', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene']], 
        chains = ['beta'], 
        organism = 'human', 
        compute_distances = False)


def test_import_vdjtools_beta_no_validation():
    import pandas as pd
    import numpy as np
    import os
    from tcrdist.paths import path_to_base
    from tcrdist.vdjtools_funcs import import_vdjtools
    vdj_tools_file_beta = os.path.join(path_to_base, 'tcrdist','data','formats','vdj.M_15_CD8_beta.clonotypes.TRB.txt.gz')
    df_beta = import_vdjtools(   vdj_tools_file = vdj_tools_file_beta ,
                        chain = 'beta',
                        organism = 'human',
                        db_file = 'alphabeta_gammadelta_db.tsv',
                        validate = False)
    assert np.all(df_beta.columns == ['count', 'freq', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'cdr3_b_nucseq','valid_v', 'valid_j', 'valid_cdr3'])
    assert False in df_beta.valid_cdr3
    assert False in df_beta.valid_v
    assert False in df_beta.valid_j


def test_import_vdjtools_alpha_w_validation():
    import pandas as pd
    import numpy as np
    import os
    from tcrdist.paths import path_to_base
    from tcrdist.vdjtools_funcs import import_vdjtools
    vdj_tools_file_alpha = os.path.join(path_to_base, 'tcrdist','data','formats','vdj.M_15_CD8_alpha.clonotypes.TRA.txt.gz')
    df_alpha = import_vdjtools( vdj_tools_file = vdj_tools_file_alpha,
                        chain = 'alpha',
                        organism = 'human',
                        db_file = 'alphabeta_gammadelta_db.tsv',
                        validate = True)
    assert np.all(df_alpha.columns == ['count', 'freq', 'cdr3_a_aa', 'v_a_gene', 'j_a_gene', 'cdr3_a_nucseq','valid_v', 'valid_j', 'valid_cdr3'])

def test_import_vdjtools_alpha_no_validation():
    import pandas as pd
    import numpy as np
    import os
    from tcrdist.paths import path_to_base
    from tcrdist.vdjtools_funcs import import_vdjtools
    vdj_tools_file_alpha = os.path.join(path_to_base, 'tcrdist','data','formats','vdj.M_15_CD8_alpha.clonotypes.TRA.txt.gz')
    df_alpha = import_vdjtools( vdj_tools_file = vdj_tools_file_alpha,
                        chain = 'alpha',
                        organism = 'human',
                        db_file = 'alphabeta_gammadelta_db.tsv',
                        validate = False)
    assert np.all(df_alpha.columns == ['count', 'freq', 'cdr3_a_aa', 'v_a_gene', 'j_a_gene', 'cdr3_a_nucseq','valid_v', 'valid_j', 'valid_cdr3'])
    assert False in df_alpha.valid_cdr3
    assert False in df_alpha.valid_v
    assert False in df_alpha.valid_j

