"""
mappers
=======

The mappers module contain dictionaries and methods for taking common output and
reformatting it to match tcrdist2 input requirements.

The top of the file contains ordered dictionaries mapping keys

More complex mappings are encapsulated in methods to be called on pd.DataFrames.

"""
import ast
import pandas as pd
from collections import OrderedDict
from functools import reduce
import warnings

tcrdist_to_tcrdist2_mapping = OrderedDict([('id', 'id'),
                                           ('epitope', 'epitope'),
                                           ('subject', 'subject'),
                                           ('cdr3a'  , 'cdr3_a_aa'),
                                           ('cdr3b'  , 'cdr3_b_aa'),
                                           ('ja_gene', 'j_a_gene'),
                                           ('va_gene', 'v_a_gene'),
                                           ('jb_gene', 'j_b_gene'),
                                           ('vb_gene', 'v_b_gene')])


vdjdb_to_tcrdist2_mapping_TRA = OrderedDict([('complex.id', 'complex_id'),
                                        ('Gene', 'gene'),
                                        ('CDR3', 'cdr3_a_aa'),
                                        ('V', 'v_a_gene'),
                                        ('J', 'j_a_gene'),
                                        ('Species', 'organism'),
                                        ('MHC A', 'mhc_a_a'),
                                        ('MHC B', 'mhc_a_b'),
                                        ('MHC class', 'mhc_a_class'),
                                        ('Epitope', 'epitope'),
                                        ('Epitope gene', 'epitope_gene'),
                                        ('Epitope species', 'epitope_species'),
                                        ('Reference', 'reference'),
                                        ('Method', 'method'),
                                        ('Meta', 'meta'),
                                        ('CDR3fix', 'cdr3fix'),
                                        ('Score', 'score')])

vdjdb_to_tcrdist2_mapping_TRB = OrderedDict([('complex.id', 'complex_id'),
                                        ('Gene', 'gene'),
                                        ('CDR3', 'cdr3_b_aa'),
                                        ('V', 'v_b_gene'),
                                        ('J', 'j_b_gene'),
                                        ('Species', 'organism'),
                                        ('MHC A', 'mhc_b_a'),
                                        ('MHC B', 'mhc_b_b'),
                                        ('MHC class', 'mhc_b_class'),
                                        ('Epitope', 'epitope'),
                                        ('Epitope gene', 'epitope_gene'),
                                        ('Epitope species', 'epitope_species'),
                                        ('Reference', 'reference'),
                                        ('Method', 'method'),
                                        ('Meta', 'meta'),
                                        ('CDR3fix', 'cdr3fix'),
                                        ('Score', 'score')])


tcrdist_clone_df_to_tcrdist2_mapping = OrderedDict([ ('clone_id'  , 'clone_id' ),
                                                     ('subject'   , 'subject'  ),
                                                     ('cdr3a'     , 'cdr3_a_aa'),
                                                     ('cdr3b'     , 'cdr3_b_aa'),
                                                     ('clone_size', 'count'    ),
                                                     ('epitope'   , 'epitope'  ),
                                                     ('ja_rep'    , 'j_a_gene' ),
                                                     ('jb_rep'    , 'j_b_gene' ),
                                                     ('subject'   , 'subject'  ),
                                                     ('va_rep'    , 'v_a_gene' ),
                                                     ('vb_rep'    , 'v_b_gene' ),
                                                     ('cdr3a_nucseq', 'cdr3_a_nucseq'),
                                                     ('cdr3b_nucseq', 'cdr3_b_nucseq'),
                                                     ('va_countreps', 'va_countreps'),
                                                     ('ja_countreps', 'ja_countreps'),
                                                     ('vb_countreps', 'vb_countreps'),
                                                     ('jb_countreps', 'jb_countreps'),
                                                     ('va_gene'     , 'va_gene'),
                                                     ('ja_gene'     , 'ja_gene'),
                                                     ('vb_gene'     , 'vb_gene'),
                                                     ('jb_gene'     , 'jb_gene')])

tcrdist2_to_tcrdist_clone_df_mapping_gd_red = OrderedDict([('clone_id', 'clone_id'),
                                                        ('subject', 'subject'),
                                                        ('cdr3_g_aa', 'cdr3a'),
                                                        ('cdr3_d_aa', 'cdr3b'),
                                                        ('count', 'clone_size'),
                                                        ('epitope', 'epitope'),
                                                        ('j_g_gene', 'ja_rep'),
                                                        ('j_d_gene', 'jb_rep'),
                                                        ('v_g_gene', 'va_rep'),
                                                        ('v_d_gene', 'vb_rep'),
                                                        ('cdr3_g_nucseq', 'cdr3a_nucseq'),
                                                        ('cdr3_d_nucseq', 'cdr3b_nucseq'),
                                                        ('vg_countreps', 'va_countreps'),
                                                        ('jg_countreps', 'ja_countreps'),
                                                        ('vd_countreps', 'vb_countreps'),
                                                        ('jd_countreps', 'jb_countreps'),
                                                        ('vg_gene'     , 'va_gene'),
                                                        ('jg_gene'     , 'ja_gene'),
                                                        ('vd_gene'     , 'vb_gene'),
                                                        ('jd_gene'     , 'jb_gene'),
                                                        ('va_countreps', 'va_countreps'),
                                                        ('ja_countreps', 'ja_countreps'),
                                                        ('vb_countreps', 'vb_countreps'),
                                                        ('jb_countreps', 'jb_countreps'),
                                                        ('va_gene'     , 'va_gene'),
                                                        ('ja_gene'     , 'ja_gene'),
                                                        ('vb_gene'     , 'vb_gene'),
                                                        ('jb_gene'     , 'jb_gene')])

tcrdist2_to_tcrdist_clone_df_mapping_gd = OrderedDict([('clone_id', 'clone_id'),
                                                        ('subject', 'subject'),
                                                        ('cdr3_g_aa', 'cdr3a'),
                                                        ('cdr3_d_aa', 'cdr3b'),
                                                        ('count', 'clone_size'),
                                                        ('epitope', 'epitope'),
                                                        ('j_g_gene', 'ja_rep'),
                                                        ('j_d_gene', 'jb_rep'),
                                                        ('v_g_gene', 'va_rep'),
                                                        ('v_d_gene', 'vb_rep'),
                                                        ('cdr3_g_nucseq', 'cdr3a_nucseq'),
                                                        ('cdr3_d_nucseq', 'cdr3b_nucseq'),
                                                        ('va_countreps', 'va_countreps'),
                                                        ('ja_countreps', 'ja_countreps'),
                                                        ('vb_countreps', 'vb_countreps'),
                                                        ('jb_countreps', 'jb_countreps'),
                                                        ('va_gene'     , 'va_gene'),
                                                        ('ja_gene'     , 'ja_gene'),
                                                        ('vb_gene'     , 'vb_gene'),
                                                        ('jb_gene'     , 'jb_gene')])



tcrdist2_to_tcrdist_clone_df_mapping = OrderedDict([('clone_id', 'clone_id'),
                                                     ('subject', 'subject'),
                                                     ('cdr3_a_aa', 'cdr3a'),
                                                     ('cdr3_b_aa', 'cdr3b'),
                                                     ('count', 'clone_size'),
                                                     ('epitope', 'epitope'),
                                                     ('j_a_gene', 'ja_rep'),
                                                     ('j_b_gene', 'jb_rep'),
                                                     ('v_a_gene', 'va_rep'),
                                                     ('v_b_gene', 'vb_rep'),
                                                     ('cdr3_a_nucseq', 'cdr3a_nucseq'),
                                                     ('cdr3_b_nucseq', 'cdr3b_nucseq'),
                                                     ('va_countreps', 'va_countreps'),
                                                     ('ja_countreps', 'ja_countreps'),
                                                     ('vb_countreps', 'vb_countreps'),
                                                     ('jb_countreps', 'jb_countreps'),
                                                     ('va_gene'     , 'va_gene'),
                                                     ('ja_gene'     , 'ja_gene'),
                                                     ('vb_gene'     , 'vb_gene'),
                                                     ('jb_gene'     , 'jb_gene')])






tcrdist2_to_tcrdist_clone_df_mapping_a = OrderedDict([('clone_id', 'clone_id'),
                                                     ('subject', 'subject'),
                                                     ('cdr3_a_aa', 'cdr3a'),
                                                     ('count', 'clone_size'),
                                                     ('epitope', 'epitope'),
                                                     ('j_a_gene', 'ja_rep'),
                                                     ('v_a_gene', 'va_rep'),
                                                     ('cdr3_a_nucseq', 'cdr3a_nucseq'),
                                                     ('cdr3_b_nucseq', 'cdr3b_nucseq'),
                                                     ('va_countreps', 'va_countreps'),
                                                     ('ja_countreps', 'ja_countreps'),
                                                     ('va_gene'     , 'va_gene'),
                                                     ('ja_gene'     , 'ja_gene')])


tcrdist2_to_tcrdist_clone_df_mapping_b = OrderedDict([('clone_id', 'clone_id'),
                                                     ('subject', 'subject'),
                                                     ('cdr3_b_aa', 'cdr3b'),
                                                     ('count', 'clone_size'),
                                                     ('epitope', 'epitope'),
                                                     ('j_b_gene', 'jb_rep'),
                                                     ('v_b_gene', 'vb_rep'),
                                                     ('cdr3_b_nucseq', 'cdr3b_nucseq'),
                                                     ('vb_countreps', 'vb_countreps'),
                                                     ('jb_countreps', 'jb_countreps'),
                                                     ('vb_gene'     , 'vb_gene'),
                                                     ('jb_gene'     , 'jb_gene')])



# Note that original tcrdist had no concept of gamma and delta syntax g->a , d->b 


tcrdist2_to_tcrdist_clone_df_mapping_g = OrderedDict([('clone_id', 'clone_id'),
                                                        ('subject', 'subject'),
                                                        ('cdr3_g_aa', 'cdr3a'),
                                                        ('count', 'clone_size'),
                                                        ('epitope', 'epitope'),
                                                        ('j_g_gene', 'ja_rep'),
                                                        ('v_g_gene', 'va_rep'),
                                                        ('cdr3_g_nucseq', 'cdr3a_nucseq'),
                                                        ('va_countreps', 'va_countreps'),
                                                        ('ja_countreps', 'ja_countreps'),
                                                        ('va_gene'     , 'va_gene'),
                                                        ('ja_gene'     , 'ja_gene')])

tcrdist2_to_tcrdist_clone_df_mapping_d = OrderedDict([('clone_id', 'clone_id'),
                                                        ('subject', 'subject'),
                                                        ('cdr3_d_aa', 'cdr3b'),
                                                        ('count', 'clone_size'),
                                                        ('epitope', 'epitope'),
                                                        ('j_d_gene', 'jb_rep'),
                                                        ('v_d_gene', 'vb_rep'),
                                                        ('cdr3_d_nucseq', 'cdr3b_nucseq'),
                                                        ('vb_countreps', 'vb_countreps'),
                                                        ('jb_countreps', 'jb_countreps'),
                                                        ('vb_gene'     , 'vb_gene'),
                                                        ('jb_gene'     , 'jb_gene')])

TCRrep_clone_df_to_TCRMotif_clone_df =OrderedDict([ ('subject'   , 'subject'),
                                                    ('epitope'  , 'epitope'),
                                                    ('v_a_gene' , 'va_rep'),
                                                    ('j_a_gene' , 'ja_rep'),
                                                    ('v_b_gene' , 'vb_rep'),
                                                    ('j_b_gene' , 'jb_rep'),
                                                    ('cdr3_a_aa', 'cdr3a'),
                                                    ('cdr3_b_aa', 'cdr3b' ) ])

TCRsubset_clone_df_to_TCRMotif_clone_df  = OrderedDict([('subject', 'subject'),
                                                        ('epitope', 'epitope'),
                                                        ('v_a_gene', 'va_rep'),
                                                        ('j_a_gene', 'ja_rep'),
                                                        ('v_b_gene', 'vb_rep'),
                                                        ('j_b_gene', 'jb_rep'),
                                                        ('cdr3_a_aa', 'cdr3a'),
                                                        ('cdr3_b_aa', 'cdr3b')])
"""gd_redundand stands for gamma delta redundant mapper"""
TCRsubset_clone_df_to_TCRMotif_clone_df_gd_red  = OrderedDict([('subject', 'subject'),
                                                                ('epitope', 'epitope'),
                                                                ('v_a_gene', 'va_rep'),
                                                                ('v_g_gene', 'va_rep'),
                                                                ('j_a_gene', 'ja_rep'),
                                                                ('j_g_gene', 'ja_rep'),
                                                                ('v_b_gene', 'vb_rep'),
                                                                ('v_d_gene', 'vb_rep'),
                                                                ('j_b_gene', 'jb_rep'),
                                                                ('j_d_gene', 'jb_rep'),
                                                                ('cdr3_a_aa', 'cdr3a'),
                                                                ('cdr3_g_aa', 'cdr3a'),
                                                                ('cdr3_b_aa', 'cdr3b'),
                                                                ('cdr3_d_aa', 'cdr3b')])


def generic_pandas_mapper(df,mapping, allow_missing = False):
    """
    Parameters
    ----------
    df : DataFrame
        Input DataFrame
    mapping : OrderedDict
        defines columns to select and rename

    Returns
    -------
    df2 : DataFrame
        Renamed DataFrame

    Example
    -------
    >>> import pandas as pd
    >>> from collections import OrderedDict
    >>> df = pd.DataFrame({"a":[1,2],"b":[3,4]})
    >>> mapper = OrderedDict([('a', 'ardvark')])
    >>> generic_pandas_mapping(df, mapper)
           ardvark
    0        1
    1        2
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a Pandas DataFrame")
    if not isinstance(mapping, OrderedDict):
        raise TypeError('mapping must be and OrderedDict')
    if not all([elem in df.keys() for elem in list(mapping.keys())]):
        if allow_missing:
            mapping = {k:v for k,v in mapping.items() if k in df.keys()}
        else:
            raise KeyError('all mapping keys must be columns of input df')

    df2 = df[mapping.keys()].rename(columns = mapping).copy()
    return(df2)

def vdjdb_to_tcrdist2(pd_df):
    """
    reformat .tsv downloaded from vdjdb (July 2019) to input format of tcrdist2 input format

    Parameters
    ----------
    pd_df : pandas.core.DataFrame likley derived from pd.read_csv() on .tsv from
    vdjdb.org

    Returns
    -------
    tcrdist2_formatted_pd_df

    Notes
    -----
    Code is written for clarity rather than length or computatinal efficiency
    since parsing demands may change.

    Note that in the unpack phase in the fist block of code a dictionary is made
    to hold all features in the vdjdb output. Much of the information from
    paired sequences is redundant (i.e. epitope). If memory demands increase
    for larger number of rows, consider modifying to unpack less of the output.

    Requires: pandas and ast (abstract syntax trees)

    Raises
    ------
    AssertionError if 'Gene' field is not TRA and TRB

    TODO
    ----
    Handle case of missing alpha or beta seq


    """
    d = {}
    # unpack vdjdb panda.DataFrame to dictionary
    for index, row in pd_df.iterrows():

        # gene will be either TRA or
        gene = row['Gene']
        assert(gene in ["TRA", "TRB"]), "Unexpected Input in vdjdb 'Gene' field not TRA or TRB. Check input."

        # complex_id is shared for paired chains
        complex_id = row['complex.id']

        # note: ast (abstract syntax trees) allows converstion of string to dictionary)
        # packed within a pd.DF cell
        meta_dict =  ast.literal_eval(row['Meta'])
        method_dict = ast.literal_eval(row['Method'])

        if gene == "TRA":
            row = row.rename(vdjdb_to_tcrdist2_mapping_TRA)
        elif gene == "TRB":
            row = row.rename(vdjdb_to_tcrdist2_mapping_TRB)

        d.setdefault(complex_id, {})[gene] = row.to_dict()
        d[complex_id][gene].update(meta_dict)
        d[complex_id][gene].update(method_dict)


    # output select fields to a list of dictionaries (l_out)
    complex_ids = sorted(d.keys())
    l_out = []
    for complex_id in complex_ids:
        try:
            id         = d[complex_id]["TRA"]['complex_id']
            cell_type  = d[complex_id]["TRA"]['cell.subset']
            organism   = d[complex_id]["TRA"]['organism']
            epitope_aa = d[complex_id]["TRA"]['epitope']
            epitope    = d[complex_id]["TRA"]['epitope.id']
            subject    = d[complex_id]["TRA"]['subject.id']

            mhc_a_a    = d[complex_id]["TRA"]['mhc_a_a']
            mhc_a_b    = d[complex_id]["TRA"]['mhc_a_b']

            mhc_b_a    = d[complex_id]["TRB"]['mhc_b_a']
            mhc_b_b    = d[complex_id]["TRB"]['mhc_b_b']

            cdr3_a_aa  = d[complex_id]["TRA"]['cdr3_a_aa']
            v_a_gene   = d[complex_id]["TRA"]['v_a_gene']
            j_a_gene   = d[complex_id]["TRA"]['j_a_gene']

            cdr3_b_aa  = d[complex_id]["TRB"]['cdr3_b_aa']
            v_b_gene   = d[complex_id]["TRB"]['v_b_gene']
            j_b_gene   = d[complex_id]["TRB"]['j_b_gene']

            frequency  = d[complex_id]["TRA"]["frequency"]
            try:
                count      = int(frequency.split("/")[0])
            except ValueError:
                count      = 1
            l_out.append(
            {'id'        : id,
            'cell_type'  : cell_type,
            'organism'   : organism,
            'epitope_aa' : epitope_aa,
            'epitope'    : epitope,
            'subject'    : subject,
            'mhc_a_a'    : mhc_a_a,
            'mhc_a_b'    : mhc_a_b,
            'mhc_b_a'    : mhc_b_a,
            'mhc_b_b'    : mhc_b_b,
            'cdr3_a_aa'  : cdr3_a_aa,
            'v_a_gene'   : v_a_gene,
            'j_a_gene'   : j_a_gene,
            'cdr3_b_aa'  : cdr3_b_aa,
            'v_b_gene'   : v_b_gene,
            'j_b_gene'   : j_b_gene,
            'frequency'  : frequency,
            'count'      : count})
        except KeyError:
            pass


    # convert list of dictionaries to pandas DataFrame
    tcrdist2_formatted_pd_df = pd.DataFrame.from_dict(l_out)
    return(tcrdist2_formatted_pd_df)



def multistudy_vdjdb_to_tcrdist2(pd_df):
    """

    map vdjdb data from multiple studies to one tcrdist2 formatted pd.dataframe

    Parameters
    ----------
    pd_df : Pandas DataFrame
        from vdjdb.org with headers as of July 10, 2019

    Returns
    -------

    tcr_formatted_dfs_merged : Pandas DataFrame
        with all paired data from input, single combined data frame tcrdist2 mapped

    """
    # 2 complex.id which links paired sequences is Reference (study specific)
    studies = pd_df['Reference'].unique()

    # break full df into sub dfs: one per study
    dfs_split_by_study = {study: pd.DataFrame for study in studies}
    for study in list(dfs_split_by_study.keys()):
        dfs_split_by_study[study] = pd_df[:][pd_df.Reference == study].reset_index()

    # tcrdist format each study specific dataframe
    tcr_formatted_dfs_split_by_study = {study: None for study in studies}
    for study in list(dfs_split_by_study.keys()):
        tcr_formatted_dfs_split_by_study[study] = vdjdb_to_tcrdist2(pd_df = dfs_split_by_study[study])

    # create list of studies with paired data
    studies_with_AB_data = [k for k in list(tcr_formatted_dfs_split_by_study.keys()) if \
                            tcr_formatted_dfs_split_by_study[k].shape[0] is not 0]
    studies_without_AB_data = [k for k in list(tcr_formatted_dfs_split_by_study.keys()) if \
                               tcr_formatted_dfs_split_by_study[k].shape[0] is not 0]

    # subset dictionary of dataframes to those with paired data
    ab = { k: tcr_formatted_dfs_split_by_study[k] for k in studies_with_AB_data }

    # count rows in each dataset for human interest
    ab_rows = [x[1].shape[0] for x in list(ab.items())]

    # convert dictionary of dataframes into a list and only keep the dataframe ignoring the study name.
    def append_study_name(name, df):
        df["reference"] = name
        return(df)

    ab_dfs = [append_study_name(x[0], x[1]) for x in list(ab.items())]

    # combine the list of dataframes into a single dataframe
    tcr_formatted_dfs_merged = reduce(lambda  x, y: pd.concat([x, y]), ab_dfs)

    return tcr_formatted_dfs_merged


def map_genes_to_first_alleles(genes, organism):
    """
    returns a list of the first alleles ('TRAJ11*01') for a given gene (TRAJ11)

    Parameters
    ----------
    genes : lists
        list of gene names using IMGT Nomenclature
    organism : string
        'human' or 'mouse'9

    Returns
    -------
    alleles : list

    Notes
    -----

    The mapping dictionary was produced with the following code
    from collections import OrderedDict
    pd_df = pd.read_csv('tcrdist/db/alphabeta_db.tsv', sep = '\t')
    g = pd_df.id.str.split("*")
    gene = [x[0] for x in g]
    X = pd.DataFrame({'gene' : gene, 'organism' : pd_df.organism.copy(), 'id' : pd_df.id.copy()})
    X = X.groupby(['gene', 'organism']).first().reset_index()
    d1 = {"mouse": OrderedDict(), "human" : OrderedDict()}
    for row in X.iterrows():
        d1[row[1]['organism']][row[1]['gene']] = row[1]['id']

    # Load second file
    pd_df = pd.read_csv('tcrdist/db/gammadelta_db.tsv', sep = '\t')
    g = pd_df.id.str.split("*")
    gene = [x[0] for x in g]
    X = pd.DataFrame({'gene' : gene, 'organism' : pd_df.organism.copy(), 'id' : pd_df.id.copy()})
    X = X.groupby(['gene', 'organism']).first().reset_index()
    for row in X.iterrows():
        d1[row[1]['organism']][row[1]['gene']] = row[1]['id']

    """

    d = {'mouse': OrderedDict([('TRAJ11', 'TRAJ11*01'), ('TRAJ12', 'TRAJ12*01'), ('TRAJ13', 'TRAJ13*01'), ('TRAJ15', 'TRAJ15*01'), ('TRAJ16', 'TRAJ16*01'), ('TRAJ17', 'TRAJ17*01'), ('TRAJ18', 'TRAJ18*01'), ('TRAJ19', 'TRAJ19*01'), ('TRAJ2', 'TRAJ2*01'), ('TRAJ20', 'TRAJ20*01'), ('TRAJ21', 'TRAJ21*01'), ('TRAJ22', 'TRAJ22*01'), ('TRAJ23', 'TRAJ23*01'), ('TRAJ24', 'TRAJ24*01'), ('TRAJ25', 'TRAJ25*01'), ('TRAJ26', 'TRAJ26*01'), ('TRAJ27', 'TRAJ27*01'), ('TRAJ28', 'TRAJ28*01'), ('TRAJ29', 'TRAJ29*01'), ('TRAJ3', 'TRAJ3*01'), ('TRAJ30', 'TRAJ30*01'), ('TRAJ31', 'TRAJ31*01'), ('TRAJ32', 'TRAJ32*01'), ('TRAJ33', 'TRAJ33*01'), ('TRAJ34', 'TRAJ34*01'), ('TRAJ35', 'TRAJ35*01'), ('TRAJ36', 'TRAJ36*01'), ('TRAJ37', 'TRAJ37*01'), ('TRAJ38', 'TRAJ38*01'), ('TRAJ39', 'TRAJ39*01'), ('TRAJ4', 'TRAJ4*01'), ('TRAJ40', 'TRAJ40*01'), ('TRAJ41', 'TRAJ41*01'), ('TRAJ42', 'TRAJ42*01'), ('TRAJ43', 'TRAJ43*01'), ('TRAJ44', 'TRAJ44*01'), ('TRAJ45', 'TRAJ45*01'), ('TRAJ46', 'TRAJ46*01'), ('TRAJ47', 'TRAJ47*01'), ('TRAJ48', 'TRAJ48*01'), ('TRAJ49', 'TRAJ49*01'), ('TRAJ5', 'TRAJ5*01'), ('TRAJ50', 'TRAJ50*01'), ('TRAJ52', 'TRAJ52*01'), ('TRAJ53', 'TRAJ53*01'), ('TRAJ54', 'TRAJ54*01'), ('TRAJ56', 'TRAJ56*01'), ('TRAJ57', 'TRAJ57*01'), ('TRAJ58', 'TRAJ58*01'), ('TRAJ59', 'TRAJ59*01'), ('TRAJ6', 'TRAJ6*01'), ('TRAJ60', 'TRAJ60*01'), ('TRAJ61', 'TRAJ61*01'), ('TRAJ7', 'TRAJ7*01'), ('TRAJ9', 'TRAJ9*01'), ('TRAV1', 'TRAV1*01'), ('TRAV10', 'TRAV10*01'), ('TRAV10D', 'TRAV10D*01'), ('TRAV10N', 'TRAV10N*01'), ('TRAV11', 'TRAV11*01'), ('TRAV11D', 'TRAV11D*01'), ('TRAV11N', 'TRAV11N*01'), ('TRAV12-1', 'TRAV12-1*01'), ('TRAV12-2', 'TRAV12-2*01'), ('TRAV12-3', 'TRAV12-3*01'), ('TRAV12D-1', 'TRAV12D-1*01'), ('TRAV12D-2', 'TRAV12D-2*01'), ('TRAV12D-3', 'TRAV12D-3*01'), ('TRAV12N-1', 'TRAV12N-1*01'), ('TRAV12N-2', 'TRAV12N-2*01'), ('TRAV12N-3', 'TRAV12N-3*01'), ('TRAV13-1', 'TRAV13-1*01'), ('TRAV13-2', 'TRAV13-2*01'), ('TRAV13-3', 'TRAV13-3*01'), ('TRAV13-4/DV7', 'TRAV13-4/DV7*01'), ('TRAV13-5', 'TRAV13-5*01'), ('TRAV13D-1', 'TRAV13D-1*01'), ('TRAV13D-2', 'TRAV13D-2*01'), ('TRAV13D-3', 'TRAV13D-3*01'), ('TRAV13D-4', 'TRAV13D-4*01'), ('TRAV13N-1', 'TRAV13N-1*01'), ('TRAV13N-2', 'TRAV13N-2*01'), ('TRAV13N-3', 'TRAV13N-3*01'), ('TRAV13N-4', 'TRAV13N-4*01'), ('TRAV14-1', 'TRAV14-1*01'), ('TRAV14-2', 'TRAV14-2*01'), ('TRAV14-3', 'TRAV14-3*01'), ('TRAV14D-1', 'TRAV14D-1*01'), ('TRAV14D-2', 'TRAV14D-2*01'), ('TRAV14D-3/DV8', 'TRAV14D-3/DV8*01'), ('TRAV14N-1', 'TRAV14N-1*01'), ('TRAV14N-2', 'TRAV14N-2*01'), ('TRAV14N-3', 'TRAV14N-3*01'), ('TRAV15-1/DV6-1', 'TRAV15-1/DV6-1*01'), ('TRAV15-2/DV6-2', 'TRAV15-2/DV6-2*01'), ('TRAV15D-1/DV6D-1', 'TRAV15D-1/DV6D-1*01'), ('TRAV15D-2/DV6D-2', 'TRAV15D-2/DV6D-2*01'), ('TRAV15N-1', 'TRAV15N-1*01'), ('TRAV15N-2', 'TRAV15N-2*01'), ('TRAV16', 'TRAV16*01'), ('TRAV16D/DV11', 'TRAV16D/DV11*01'), ('TRAV16N', 'TRAV16N*01'), ('TRAV17', 'TRAV17*01'), ('TRAV18', 'TRAV18*01'), ('TRAV19', 'TRAV19*01'), ('TRAV2', 'TRAV2*01'), ('TRAV20', 'TRAV20*01'), ('TRAV21/DV12', 'TRAV21/DV12*01'), ('TRAV3-1', 'TRAV3-1*01'), ('TRAV3-3', 'TRAV3-3*01'), ('TRAV3-4', 'TRAV3-4*01'), ('TRAV3D-3', 'TRAV3D-3*01'), ('TRAV3N-3', 'TRAV3N-3*01'), ('TRAV4-2', 'TRAV4-2*01'), ('TRAV4-3', 'TRAV4-3*01'), ('TRAV4-4/DV10', 'TRAV4-4/DV10*01'), ('TRAV4D-2', 'TRAV4D-2*01'), ('TRAV4D-3', 'TRAV4D-3*01'), ('TRAV4D-4', 'TRAV4D-4*01'), ('TRAV4N-3', 'TRAV4N-3*01'), ('TRAV4N-4', 'TRAV4N-4*01'), ('TRAV5-1', 'TRAV5-1*01'), ('TRAV5-2', 'TRAV5-2*01'), ('TRAV5-4', 'TRAV5-4*01'), ('TRAV5D-2', 'TRAV5D-2*01'), ('TRAV5D-4', 'TRAV5D-4*01'), ('TRAV5N-2', 'TRAV5N-2*01'), ('TRAV5N-4', 'TRAV5N-4*01'), ('TRAV6-1', 'TRAV6-1*01'), ('TRAV6-2', 'TRAV6-2*01'), ('TRAV6-3', 'TRAV6-3*01'), ('TRAV6-4', 'TRAV6-4*01'), ('TRAV6-5', 'TRAV6-5*01'), ('TRAV6-6', 'TRAV6-6*01'), ('TRAV6-7/DV9', 'TRAV6-7/DV9*01'), ('TRAV6D-3', 'TRAV6D-3*01'), ('TRAV6D-4', 'TRAV6D-4*01'), ('TRAV6D-5', 'TRAV6D-5*01'), ('TRAV6D-6', 'TRAV6D-6*01'), ('TRAV6D-7', 'TRAV6D-7*01'), ('TRAV6N-5', 'TRAV6N-5*01'), ('TRAV6N-6', 'TRAV6N-6*01'), ('TRAV6N-7', 'TRAV6N-7*01'), ('TRAV7-1', 'TRAV7-1*01'), ('TRAV7-2', 'TRAV7-2*01'), ('TRAV7-3', 'TRAV7-3*01'), ('TRAV7-4', 'TRAV7-4*01'), ('TRAV7-5', 'TRAV7-5*01'), ('TRAV7-6', 'TRAV7-6*01'), ('TRAV7D-2', 'TRAV7D-2*01'), ('TRAV7D-3', 'TRAV7D-3*01'), ('TRAV7D-4', 'TRAV7D-4*01'), ('TRAV7D-5', 'TRAV7D-5*01'), ('TRAV7D-6', 'TRAV7D-6*01'), ('TRAV7N-4', 'TRAV7N-4*01'), ('TRAV7N-5', 'TRAV7N-5*01'), ('TRAV7N-6', 'TRAV7N-6*01'), ('TRAV8-1', 'TRAV8-1*01'), ('TRAV8-2', 'TRAV8-2*01'), ('TRAV8D-1', 'TRAV8D-1*01'), ('TRAV8D-2', 'TRAV8D-2*01'), ('TRAV8N-2', 'TRAV8N-2*01'), ('TRAV9-1', 'TRAV9-1*01'), ('TRAV9-2', 'TRAV9-2*01'), ('TRAV9-3', 'TRAV9-3*01'), ('TRAV9-4', 'TRAV9-4*01'), ('TRAV9D-1', 'TRAV9D-1*01'), ('TRAV9D-2', 'TRAV9D-2*01'), ('TRAV9D-3', 'TRAV9D-3*01'), ('TRAV9D-4', 'TRAV9D-4*01'), ('TRAV9N-2', 'TRAV9N-2*01'), ('TRAV9N-3', 'TRAV9N-3*01'), ('TRAV9N-4', 'TRAV9N-4*01'), ('TRBD1', 'TRBD1*01'), ('TRBD2', 'TRBD2*01'), ('TRBJ1-1', 'TRBJ1-1*01'), ('TRBJ1-2', 'TRBJ1-2*01'), ('TRBJ1-3', 'TRBJ1-3*01'), ('TRBJ1-4', 'TRBJ1-4*01'), ('TRBJ1-5', 'TRBJ1-5*01'), ('TRBJ1-6', 'TRBJ1-6*01'), ('TRBJ1-7', 'TRBJ1-7*01'), ('TRBJ2-1', 'TRBJ2-1*01'), ('TRBJ2-2', 'TRBJ2-2*01'), ('TRBJ2-3', 'TRBJ2-3*01'), ('TRBJ2-4', 'TRBJ2-4*01'), ('TRBJ2-5', 'TRBJ2-5*01'), ('TRBJ2-6', 'TRBJ2-6*01'), ('TRBJ2-7', 'TRBJ2-7*01'), ('TRBV1', 'TRBV1*01'), ('TRBV10', 'TRBV10*01'), ('TRBV12-1', 'TRBV12-1*01'), ('TRBV12-2', 'TRBV12-2*01'), ('TRBV13-1', 'TRBV13-1*01'), ('TRBV13-2', 'TRBV13-2*01'), ('TRBV13-3', 'TRBV13-3*01'), ('TRBV14', 'TRBV14*01'), ('TRBV15', 'TRBV15*01'), ('TRBV16', 'TRBV16*01'), ('TRBV17', 'TRBV17*01'), ('TRBV19', 'TRBV19*01'), ('TRBV2', 'TRBV2*01'), ('TRBV20', 'TRBV20*01'), ('TRBV21', 'TRBV21*01'), ('TRBV23', 'TRBV23*01'), ('TRBV24', 'TRBV24*01'), ('TRBV26', 'TRBV26*01'), ('TRBV29', 'TRBV29*01'), ('TRBV3', 'TRBV3*01'), ('TRBV30', 'TRBV30*01'), ('TRBV31', 'TRBV31*01'), ('TRBV4', 'TRBV4*01'), ('TRBV5', 'TRBV5*01'), ('TRBV8', 'TRBV8*01'), ('TRBV9', 'TRBV9*01'), ('TRDD1', 'TRDD1*01'), ('TRDD2', 'TRDD2*01'), ('TRDJ1', 'TRDJ1*01'), ('TRDJ2', 'TRDJ2*01'), ('TRDV1', 'TRDV1*01'), ('TRDV2-1', 'TRDV2-1*01'), ('TRDV2-2', 'TRDV2-2*01'), ('TRDV4', 'TRDV4*01'), ('TRDV5', 'TRDV5*01'), ('TRGJ1', 'TRGJ1*01'), ('TRGJ2', 'TRGJ2*01'), ('TRGJ3', 'TRGJ3*01'), ('TRGJ4', 'TRGJ4*01'), ('TRGV1', 'TRGV1*01'), ('TRGV2', 'TRGV2*01'), ('TRGV3', 'TRGV3*01'), ('TRGV4', 'TRGV4*01'), ('TRGV5', 'TRGV5*01'), ('TRGV6', 'TRGV6*01'), ('TRGV7', 'TRGV7*01')]), 'human': OrderedDict([('TRAJ1', 'TRAJ1*01'), ('TRAJ10', 'TRAJ10*01'), ('TRAJ11', 'TRAJ11*01'), ('TRAJ12', 'TRAJ12*01'), ('TRAJ13', 'TRAJ13*01'), ('TRAJ14', 'TRAJ14*01'), ('TRAJ15', 'TRAJ15*01'), ('TRAJ16', 'TRAJ16*01'), ('TRAJ17', 'TRAJ17*01'), ('TRAJ18', 'TRAJ18*01'), ('TRAJ19', 'TRAJ19*01'), ('TRAJ2', 'TRAJ2*01'), ('TRAJ20', 'TRAJ20*01'), ('TRAJ21', 'TRAJ21*01'), ('TRAJ22', 'TRAJ22*01'), ('TRAJ23', 'TRAJ23*01'), ('TRAJ24', 'TRAJ24*01'), ('TRAJ25', 'TRAJ25*01'), ('TRAJ26', 'TRAJ26*01'), ('TRAJ27', 'TRAJ27*01'), ('TRAJ28', 'TRAJ28*01'), ('TRAJ29', 'TRAJ29*01'), ('TRAJ3', 'TRAJ3*01'), ('TRAJ30', 'TRAJ30*01'), ('TRAJ31', 'TRAJ31*01'), ('TRAJ32', 'TRAJ32*01'), ('TRAJ33', 'TRAJ33*01'), ('TRAJ34', 'TRAJ34*01'), ('TRAJ35', 'TRAJ35*01'), ('TRAJ36', 'TRAJ36*01'), ('TRAJ37', 'TRAJ37*01'), ('TRAJ38', 'TRAJ38*01'), ('TRAJ39', 'TRAJ39*01'), ('TRAJ4', 'TRAJ4*01'), ('TRAJ40', 'TRAJ40*01'), ('TRAJ41', 'TRAJ41*01'), ('TRAJ42', 'TRAJ42*01'), ('TRAJ43', 'TRAJ43*01'), ('TRAJ44', 'TRAJ44*01'), ('TRAJ45', 'TRAJ45*01'), ('TRAJ46', 'TRAJ46*01'), ('TRAJ47', 'TRAJ47*01'), ('TRAJ48', 'TRAJ48*01'), ('TRAJ49', 'TRAJ49*01'), ('TRAJ5', 'TRAJ5*01'), ('TRAJ50', 'TRAJ50*01'), ('TRAJ51', 'TRAJ51*01'), ('TRAJ52', 'TRAJ52*01'), ('TRAJ53', 'TRAJ53*01'), ('TRAJ54', 'TRAJ54*01'), ('TRAJ55', 'TRAJ55*01'), ('TRAJ56', 'TRAJ56*01'), ('TRAJ57', 'TRAJ57*01'), ('TRAJ58', 'TRAJ58*01'), ('TRAJ59', 'TRAJ59*01'), ('TRAJ6', 'TRAJ6*01'), ('TRAJ60', 'TRAJ60*01'), ('TRAJ61', 'TRAJ61*01'), ('TRAJ7', 'TRAJ7*01'), ('TRAJ8', 'TRAJ8*01'), ('TRAJ9', 'TRAJ9*01'), ('TRAV1-1', 'TRAV1-1*01'), ('TRAV1-2', 'TRAV1-2*01'), ('TRAV10', 'TRAV10*01'), ('TRAV11', 'TRAV11*01'), ('TRAV12-1', 'TRAV12-1*01'), ('TRAV12-2', 'TRAV12-2*01'), ('TRAV12-3', 'TRAV12-3*01'), ('TRAV13-1', 'TRAV13-1*01'), ('TRAV13-2', 'TRAV13-2*01'), ('TRAV14/DV4', 'TRAV14/DV4*01'), ('TRAV16', 'TRAV16*01'), ('TRAV17', 'TRAV17*01'), ('TRAV18', 'TRAV18*01'), ('TRAV19', 'TRAV19*01'), ('TRAV2', 'TRAV2*01'), ('TRAV20', 'TRAV20*01'), ('TRAV21', 'TRAV21*01'), ('TRAV22', 'TRAV22*01'), ('TRAV23/DV6', 'TRAV23/DV6*01'), ('TRAV24', 'TRAV24*01'), ('TRAV25', 'TRAV25*01'), ('TRAV26-1', 'TRAV26-1*01'), ('TRAV26-2', 'TRAV26-2*01'), ('TRAV27', 'TRAV27*01'), ('TRAV29/DV5', 'TRAV29/DV5*01'), ('TRAV3', 'TRAV3*01'), ('TRAV30', 'TRAV30*01'), ('TRAV34', 'TRAV34*01'), ('TRAV35', 'TRAV35*01'), ('TRAV36/DV7', 'TRAV36/DV7*01'), ('TRAV38-1', 'TRAV38-1*01'), ('TRAV38-2/DV8', 'TRAV38-2/DV8*01'), ('TRAV39', 'TRAV39*01'), ('TRAV4', 'TRAV4*01'), ('TRAV40', 'TRAV40*01'), ('TRAV41', 'TRAV41*01'), ('TRAV5', 'TRAV5*01'), ('TRAV6', 'TRAV6*01'), ('TRAV7', 'TRAV7*01'), ('TRAV8-1', 'TRAV8-1*01'), ('TRAV8-2', 'TRAV8-2*01'), ('TRAV8-3', 'TRAV8-3*01'), ('TRAV8-4', 'TRAV8-4*01'), ('TRAV8-6', 'TRAV8-6*01'), ('TRAV8-7', 'TRAV8-7*01'), ('TRAV9-1', 'TRAV9-1*01'), ('TRAV9-2', 'TRAV9-2*01'), ('TRBD1', 'TRBD1*01'), ('TRBD2', 'TRBD2*01'), ('TRBJ1-1', 'TRBJ1-1*01'), ('TRBJ1-2', 'TRBJ1-2*01'), ('TRBJ1-3', 'TRBJ1-3*01'), ('TRBJ1-4', 'TRBJ1-4*01'), ('TRBJ1-5', 'TRBJ1-5*01'), ('TRBJ1-6', 'TRBJ1-6*01'), ('TRBJ2-1', 'TRBJ2-1*01'), ('TRBJ2-2', 'TRBJ2-2*01'), ('TRBJ2-2P', 'TRBJ2-2P*01'), ('TRBJ2-3', 'TRBJ2-3*01'), ('TRBJ2-4', 'TRBJ2-4*01'), ('TRBJ2-5', 'TRBJ2-5*01'), ('TRBJ2-6', 'TRBJ2-6*01'), ('TRBJ2-7', 'TRBJ2-7*01'), ('TRBV1', 'TRBV1*01'), ('TRBV10-1', 'TRBV10-1*01'), ('TRBV10-2', 'TRBV10-2*01'), ('TRBV10-3', 'TRBV10-3*01'), ('TRBV11-1', 'TRBV11-1*01'), ('TRBV11-2', 'TRBV11-2*01'), ('TRBV11-3', 'TRBV11-3*01'), ('TRBV12-1', 'TRBV12-1*01'), ('TRBV12-2', 'TRBV12-2*01'), ('TRBV12-3', 'TRBV12-3*01'), ('TRBV12-4', 'TRBV12-4*01'), ('TRBV12-5', 'TRBV12-5*01'), ('TRBV13', 'TRBV13*01'), ('TRBV14', 'TRBV14*01'), ('TRBV15', 'TRBV15*01'), ('TRBV16', 'TRBV16*01'), ('TRBV17', 'TRBV17*01'), ('TRBV18', 'TRBV18*01'), ('TRBV19', 'TRBV19*01'), ('TRBV2', 'TRBV2*01'), ('TRBV20-1', 'TRBV20-1*01'), ('TRBV20/OR9-2', 'TRBV20/OR9-2*01'), ('TRBV21-1', 'TRBV21-1*01'), ('TRBV21/OR9-2', 'TRBV21/OR9-2*01'), ('TRBV23-1', 'TRBV23-1*01'), ('TRBV23/OR9-2', 'TRBV23/OR9-2*01'), ('TRBV24-1', 'TRBV24-1*01'), ('TRBV24/OR9-2', 'TRBV24/OR9-2*01'), ('TRBV25-1', 'TRBV25-1*01'), ('TRBV26', 'TRBV26*01'), ('TRBV26/OR9-2', 'TRBV26/OR9-2*01'), ('TRBV27', 'TRBV27*01'), ('TRBV28', 'TRBV28*01'), ('TRBV29-1', 'TRBV29-1*01'), ('TRBV29/OR9-2', 'TRBV29/OR9-2*01'), ('TRBV3-1', 'TRBV3-1*01'), ('TRBV3-2', 'TRBV3-2*01'), ('TRBV30', 'TRBV30*01'), ('TRBV4-1', 'TRBV4-1*01'), ('TRBV4-2', 'TRBV4-2*01'), ('TRBV4-3', 'TRBV4-3*01'), ('TRBV5-1', 'TRBV5-1*01'), ('TRBV5-3', 'TRBV5-3*01'), ('TRBV5-4', 'TRBV5-4*01'), ('TRBV5-5', 'TRBV5-5*01'), ('TRBV5-6', 'TRBV5-6*01'), ('TRBV5-7', 'TRBV5-7*01'), ('TRBV5-8', 'TRBV5-8*01'), ('TRBV6-1', 'TRBV6-1*01'), ('TRBV6-2', 'TRBV6-2*01'), ('TRBV6-3', 'TRBV6-3*01'), ('TRBV6-4', 'TRBV6-4*01'), ('TRBV6-5', 'TRBV6-5*01'), ('TRBV6-6', 'TRBV6-6*01'), ('TRBV6-7', 'TRBV6-7*01'), ('TRBV6-8', 'TRBV6-8*01'), ('TRBV6-9', 'TRBV6-9*01'), ('TRBV7-1', 'TRBV7-1*01'), ('TRBV7-2', 'TRBV7-2*01'), ('TRBV7-3', 'TRBV7-3*01'), ('TRBV7-4', 'TRBV7-4*01'), ('TRBV7-6', 'TRBV7-6*01'), ('TRBV7-7', 'TRBV7-7*01'), ('TRBV7-8', 'TRBV7-8*01'), ('TRBV7-9', 'TRBV7-9*01'), ('TRBV9', 'TRBV9*01'), ('TRDD1', 'TRDD1*01'), ('TRDD2', 'TRDD2*01'), ('TRDD3', 'TRDD3*01'), ('TRDJ1', 'TRDJ1*01'), ('TRDJ2', 'TRDJ2*01'), ('TRDJ3', 'TRDJ3*01'), ('TRDJ4', 'TRDJ4*01'), ('TRDV1', 'TRDV1*01'), ('TRDV2', 'TRDV2*01'), ('TRDV3', 'TRDV3*01'), ('TRGJ1', 'TRGJ1*01'), ('TRGJ2', 'TRGJ2*01'), ('TRGJP', 'TRGJP*01'), ('TRGJP1', 'TRGJP1*01'), ('TRGJP2', 'TRGJP2*01'), ('TRGV1', 'TRGV1*01'), ('TRGV10', 'TRGV10*01'), ('TRGV11', 'TRGV11*01'), ('TRGV2', 'TRGV2*01'), ('TRGV3', 'TRGV3*01'), ('TRGV4', 'TRGV4*01'), ('TRGV5', 'TRGV5*01'), ('TRGV5P', 'TRGV5P*01'), ('TRGV8', 'TRGV8*01'), ('TRGV9', 'TRGV9*01'), ('TRGVA', 'TRGVA*01')])}

    alleles = []
    for g in genes:
        try:
            alleles.append(d[organism][g])
        except KeyError:
            alleles.append(None)
            warnings.warn("{} not found in mapping DB, no allele could be mapped for this gene".format(g))
    return(alleles)

def populate_legacy_fields(df, chains =['alpha', 'beta']):
    """
    For instances when we only have v_x_gene and v_j_gene and we need to supply
    'va_countreps' and 'va_gene' for use of tcrdist.subset.TCRsubset()
    """
    if 'alpha' in chains:
        df['va_countreps'] = df['v_a_gene'].copy()
        df['ja_countreps'] = df['j_a_gene'].copy()
        df['va_gene'] = df['v_a_gene'].copy()
        df['ja_gene'] = df['j_a_gene'].copy()
    if 'beta' in chains:
        df['vb_countreps'] = df['v_b_gene'].copy()
        df['jb_countreps'] = df['j_b_gene'].copy()
        df['vb_gene'] = df['v_b_gene'].copy()
        df['jb_gene'] = df['j_b_gene'].copy()
    if 'gamma' in chains:
        df['va_countreps'] = df['v_g_gene'].copy()
        df['ja_countreps'] = df['j_g_gene'].copy()
        df['va_gene'] = df['v_g_gene'].copy()
        df['ja_gene'] = df['j_g_gene'].copy()
    if 'delta' in chains:
        df['vb_countreps'] = df['v_d_gene'].copy()
        df['jb_countreps'] = df['j_d_gene'].copy()
        df['vb_gene'] = df['v_d_gene'].copy()
        df['jb_gene'] = df['j_d_gene'].copy()
    return(df)