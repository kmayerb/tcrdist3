import os
import re
import pandas as pd
import numpy as np
import warnings
from progress.bar import IncrementalBar
from tcrdist.pgen import OlgaModel

class TCRtree():
    def __init__(self,tcrrep, html_name):
        """
        Automatic Hierarchical Cluster Plotting.

        This is a class that automates the process of producing 
        interactive TCR hierarchical clustering tree.

        These trees can be constructed manually see `docs
        <https://tcrdist3.readthedocs.io/en/latest/motif_gallery.html>_.
        <http://www.python.org/>`_.
        Parameters
        ----------

        html_name : str
            name for html file output e.g., 'DEFAULT.html'
        pwmat_str_b : str
            name of pairwise matrix attribute to be used for clustering beta chains e.g., 'pw_beta'
        pwmat_str_a : str
            name of pairwise matrix attribute to be used for clustering alpha chains e.g., 'pw_alpha'
        single : bool
            If true, make summary based on each clone being present in single-copy, 
            otherwise, the 'count' column is used when calculating percentages.
            NOTE: 'count_col' in default_hcluster_diff_kwargs can also set to 
            'single'. If true, diversity metrics will also be based on clones 
            rather than clonal abundances. 
        generate_svgs : bool 
            If True, SVG logos are produced for each node where .hcluster_df['prune'] is 0
        verbose : bool
            report on status


        default_hcluster_diff_kwargs: dict
            kwargs dictionary for (tcrdist.rep_diff.hcluster_diff)

                clone_df : pd.DataFrame [nclones x metadata]
                    Contains metadata for each clone.
                pwmat : np.ndarray [nclones x nclones]
                    Square distance matrix for defining neighborhoods
                x_cols : list
                    List of columns to be tested for association with the neighborhood
                count_col : str
                    Column in clone_df that specifies counts.
                    Default none assumes count of 1 cell for each row.
                subset_ind : None or np.ndarray with partial index of df, optional
                    Provides option to tally counts only within a subset of df, but to maintain the clustering
                    of all individuals. Allows for one clustering of pooled TCRs,
                    but tallying/testing within a subset (e.g. participants or conditions)
                hclust_method : str
                    Method for hierarchical clustering, passed to the scipy.clustering.hierarchy
                    linkage function.
                optimal_ordering : bool
                    Flag passed to the scipy.clustering.hierarchy linkage function to improve
                    visual tree layout. Can be slow for large trees.
                test_method : str or None
                    Specifies Fisher's exact test ("fishers"), Chi-squared ("chi2") or
                    Cochran-Mantel-Haenszel test ("chm") for testing.
        default_member_summ_kwargs: dict
            kwargs dictionary for (tcrdist.rep_diff.member_summ)

                Return additional summary info about each result (row)) based on the
                members of the cluster. This is helpful for preparing strings 
                to add to the tooltip in hierdiff.plot_hclust_props.

                res_df : pd.DataFrame [nclusters x result cols]
                    Returned from neighborhood_diff or hcluster_diff
                clone_df : pd.DataFrame [nclones x metadata]
                    Contains metadata for each clone.
                key_col : str
                    Column in res_df that specifies the iloc of members in the clone_df
                count_col : str
                    Column in clone_df that specifies counts.
                    Default none assumes count of 1 cell for each row.
                addl_cols : list
                    Columns to summarize
                addl_n : int
                    Number of top N clones to include in the summary of
                    each cluster.
        default_plot_hclust_props : dict
            kwargs dictionary for (hierdiff.plot_hclust_props)

                Plot tree of linkage-based hierarchical clustering, with nodes colored using stacked bars
                representing proportion of cluster members associated with specific conditions. Nodes also optionally
                annotated with pvalue, number of members or cluster ID.
                
                Z : linkage matrix
                    Result of calling sch.linkage on a compressed pair-wise distance matrix
                res : pd.DataFrame
                    Result from calling hcluster_diff, with observed/frequencies and p-values for each node
                alpha_col : str
                    Column in res to use for 'alpha' annotation
                alpha : float
                    Threshold for plotting the stacked bars and annotation
                colors : tuple of valid colors
                    Used for stacked bars of conditions at each node
                prune_col : str/column in res
                    Column of res that indicates whether a result/node can be pruned from the tree.
                    The tree will not print any pruned nodes that only contain other pruned nodes.
        """
        self.tcrrep = tcrrep
        self.html_name = html_name

        kws1,kws2,kws3 = _get_default_kwargs(chains = self.tcrrep.chains)
        self.default_hcluster_diff_kwargs = kws1
        self.default_member_summ_kwargs   = kws2
        self.default_plot_hclust_props    = kws3
        self.prune = 3
        self.combine_olga = True
        self.pwmat_str_b = 'pw_beta'
        self.pwmat_str_a = 'pw_alpha'
        self.single = True
        self.generate_svgs = True
        self.verbose = True

    def build_tree( self ):

        print(self.default_plot_hclust_props['tooltip_cols'])
        _auto_hdiff2(   tcrrep                       = self.tcrrep,
                        html_name                    = self.html_name,
                        pwmat_str_b                  = self.pwmat_str_b,
                        pwmat_str_a                  = self.pwmat_str_a,
                        single                       = self.single,
                        generate_svgs                = self.generate_svgs ,
                        combine_olga                 = self.combine_olga,
                        prune                        = self.prune,
                        verbose                      = self.verbose,
                        default_hcluster_diff_kwargs = self.default_hcluster_diff_kwargs,
                        default_member_summ_kwargs   = self.default_member_summ_kwargs,
                        default_plot_hclust_props    = self.default_plot_hclust_props)



def _get_default_kwargs(chains = ['beta']):
    default_hcluster_diff_kwargs = {
        'clone_df' : None, 
        'pwmat' : None, 
        'x_cols' : ['epitope'],#None,#['epitope'], 
        'Z' : None, 
        'count_col' : 'count', 
        'subset_ind' : None, 
        'hclust_method' : 'complete', 
        'optimal_ordering' : True, 
        'test_method' : 'fishers'
        }
    default_member_summ_kwargs = {
        'key_col' : 'neighbors_i', 
        'count_col' : 'count', 
        'addl_cols' :['subject'], 
        'addl_n' : 1
        }

    if len(chains) == 1:
        chain = chains[0]
        chain_abr = chain[0]
        default_plot_hclust_props = {
            'title':'', 
            'alpha_col':'pvalue', 
            'alpha':0.05, 
            'tooltip_cols':['subject',
                            'mean_dist',
                            'pct_dist_75',
                            'pct_dist_50',
                            'pct_dist_25',
                            'fuzzy_simpson_diversity_75', 
                            'fuzzy_simpson_diversity_50',
                            'fuzzy_simpson_diversity_25',
                            f'cdr3_{chain_abr}_aa', 
                            f'v_{chain_abr}_gene',
                            f'j_{chain_abr}_gene',
                            f'svg_{chain}',
                            f'svg_raw_{chain}',
                            f'ref_size_{chain}', 
                            f'ref_unique_{chain}',
                            f'percent_missing_{chain}']}
    elif len(chains) == 2:
        chain = chains[0]
        chain_abr = chain[0]
        chain2 = chains[1]
        chain2_abr = chain2[0]
        default_plot_hclust_props = {
            'title':'', 
            'alpha_col':'pvalue', 
            'alpha':0.05, 
            'tooltip_cols':[
                            'mean_dist',
                            'pct_dist_75',
                            'pct_dist_50',
                            'pct_dist_25',
                            'fuzzy_simpson_diversity_75', 
                            'fuzzy_simpson_diversity_50',
                            'fuzzy_simpson_diversity_25',
                            f'cdr3_{chain_abr}_aa', 
                            f'v_{chain_abr}_gene',
                            f'j_{chain_abr}_gene',
                            f'svg_{chain}',
                            f'svg_raw_{chain}',
                            f'ref_size_{chain}', 
                            f'ref_unique_{chain}',
                            f'percent_missing_{chain}',
                            f'cdr3_{chain2_abr}_aa', 
                            f'v_{chain2_abr}_gene',
                            f'j_{chain2_abr}_gene',
                            f'svg_{chain2}',
                            f'svg_raw_{chain2}',
                            f'ref_size_{chain2}', 
                            f'ref_unique_{chain2}',
                            f'percent_missing_{chain2}']}
    else:
        raise ValueError('Only 1 or 2 chains can be provided')

    return default_hcluster_diff_kwargs, default_member_summ_kwargs, default_plot_hclust_props


def _auto_hdiff2(tcrrep,
                html_name = 'DEFAULT.html',
                pwmat_str_b = 'pw_beta',
                pwmat_str_a = 'pw_alpha',
                single = True,
                generate_svgs = True,
                combine_olga = False,
                verbose = True,
                prune = 3,
                default_hcluster_diff_kwargs = _get_default_kwargs(chains = ['beta'])[0],
                default_member_summ_kwargs   = _get_default_kwargs(chains = ['beta'])[1],
                default_plot_hclust_props    = _get_default_kwargs(chains = ['beta'])[2]):
    """
    Automatic Hierarchical Cluster Plotting

    Parameters
    ----------

    html_name : str
        name for html file output e.g., 'DEFAULT.html'
    pwmat_str_b : str
        name of pairwise matrix attribute to be used for clustering beta chains e.g., 'pw_beta'
    pwmat_str_a : str
        name of pairwise matrix attribute to be used for clustering alpha chains e.g., 'pw_alpha'
    single : bool
        If true, make summary based on each clone being present in single-copy, 
        otherwise, the 'count' column is used when calculating percentages.
        NOTE: 'count_col' in default_hcluster_diff_kwargs can also set to 
        'single'. If true, diversity metrics will also be based on clones 
        rather than clonal abundances. 
    generate_svgs : bool 
        If True, SVG logos are produced for each node where .hcluster_df['prune'] is 0
    verbose : bool
        report on status
    default_hcluster_diff_kwargs: dict
        kwargs dictionary for (tcrdist.rep_diff.hcluster_diff)

            clone_df : pd.DataFrame [nclones x metadata]
                Contains metadata for each clone.
            pwmat : np.ndarray [nclones x nclones]
                Square distance matrix for defining neighborhoods
            x_cols : list
                List of columns to be tested for association with the neighborhood
            count_col : str
                Column in clone_df that specifies counts.
                Default none assumes count of 1 cell for each row.
            subset_ind : None or np.ndarray with partial index of df, optional
                Provides option to tally counts only within a subset of df, but to maintain the clustering
                of all individuals. Allows for one clustering of pooled TCRs,
                but tallying/testing within a subset (e.g. participants or conditions)
            hclust_method : str
                Method for hierarchical clustering, passed to the scipy.clustering.hierarchy
                linkage function.
            optimal_ordering : bool
                Flag passed to the scipy.clustering.hierarchy linkage function to improve
                visual tree layout. Can be slow for large trees.
            test_method : str or None
                Specifies Fisher's exact test ("fishers"), Chi-squared ("chi2") or
                Cochran-Mantel-Haenszel test ("chm") for testing.
    default_member_summ_kwargs: dict
        kwargs dictionary for (tcrdist.rep_diff.member_summ)

            Return additional summary info about each result (row)) based on the
            members of the cluster. This is helpful for preparing strings 
            to add to the tooltip in hierdiff.plot_hclust_props.

            res_df : pd.DataFrame [nclusters x result cols]
                Returned from neighborhood_diff or hcluster_diff
            clone_df : pd.DataFrame [nclones x metadata]
                Contains metadata for each clone.
            key_col : str
                Column in res_df that specifies the iloc of members in the clone_df
            count_col : str
                Column in clone_df that specifies counts.
                Default none assumes count of 1 cell for each row.
            addl_cols : list
                Columns to summarize
            addl_n : int
                Number of top N clones to include in the summary of
                each cluster.
    default_plot_hclust_props : dict
        kwargs dictionary for (hierdiff.plot_hclust_props)

            Plot tree of linkage-based hierarchical clustering, with nodes colored using stacked bars
            representing proportion of cluster members associated with specific conditions. Nodes also optionally
            annotated with pvalue, number of members or cluster ID.
            
            Z : linkage matrix
                Result of calling sch.linkage on a compressed pair-wise distance matrix
            res : pd.DataFrame
                Result from calling hcluster_diff, with observed/frequencies and p-values for each node
            alpha_col : str
                Column in res to use for 'alpha' annotation
            alpha : float
                Threshold for plotting the stacked bars and annotation
            colors : tuple of valid colors
                Used for stacked bars of conditions at each node
            prune_col : str/column in res
                Column of res that indicates whether a result/node can be pruned from the tree.
                The tree will not print any pruned nodes that only contain other pruned nodes.
    """
    import os
    import pandas as pd
    import numpy as np
    import warnings
    from tcrsampler.sampler import TCRsampler
    from palmotif import compute_pal_motif, svg_logo
    from tcrdist.adpt_funcs import get_centroid_seq, get_centroid_seq_alpha
    from tcrdist.summarize import _select
    from tcrdist.repertoire import TCRrep
    from tcrdist.rep_diff import hcluster_diff, member_summ
    from tcrdist.summarize import _select
    from tcrdist.pgen import OlgaModel
    from palmotif import compute_pal_motif, svg_logo
    from hierdiff import plot_hclust_props
    from numpy.random import randint
    from tcrdist.diversity import generalized_simpsons_entropy
    from tcrdist.diversity import fuzzy_diversity
    

    # Load clone_df directly from input object
    if default_hcluster_diff_kwargs['clone_df'] is None:
        default_hcluster_diff_kwargs['clone_df'] = getattr(tcrrep, 'clone_df')

    # Get appropirate pwmat_str (pw_beta, pw_alpha, or one of the CDRs, e.g., pw_cdr3_b_aa)
    if 'alpha' in tcrrep.chains:
        pwmat_str = pwmat_str_a
    if 'beta' in tcrrep.chains:
        pwmat_str = pwmat_str_b


    # Get appropirate pairwise matrix
    
    if default_hcluster_diff_kwargs['pwmat'] is None:
        if verbose : print(f"pwmat set with {pwmat_str}")
        default_hcluster_diff_kwargs['pwmat'] = getattr(tcrrep, pwmat_str)
    else: 
        if verbose : print("pwmat was directly provided as a kwarg")
    """Handle the Fact that 2 or more catagoorical levels are need to run hierdiff"""
    
    x_cols = default_hcluster_diff_kwargs['x_cols']
    if x_cols is None:
        tcrrep.clone_df['dummy'] = \
            [['X1','X2'][randint(2)] for x in range(tcrrep.clone_df.shape[0])]
        default_hcluster_diff_kwargs['x_cols'] = ['dummy']
        default_hcluster_diff_kwargs['test_method'] = 'fishers' 
        warnings.warn(f"Because x_cols was None, setting random dummy values, and using {default_hcluster_diff_kwargs['test_method']}\n")
    
    elif tcrrep.clone_df[x_cols].nunique()[0] == 2: 
        default_hcluster_diff_kwargs['test_method'] = 'fishers' 
    
    elif tcrrep.clone_df[x_cols].nunique()[0] > 2: 
        default_hcluster_diff_kwargs['test_method'] = 'chi2'
    
    elif tcrrep.clone_df[x_cols].nunique()[0] < 2 : 
        tcrrep.clone_df['dummy'] = \
            [['X1','X2'][randint(2)] for x in range(tcrrep.clone_df.shape[0])]
        default_hcluster_diff_kwargs['x_cols'] = ['dummy']
        default_hcluster_diff_kwargs['test_method'] = 'fishers' 
        warnings.warn(f"Because x_cols was None, setting random dummy values, and using {default_hcluster_diff_kwargs['test_method']}\n")

    """ Run hcluster_df """
    
    bar = IncrementalBar(f'Run hcluster_diff :', max = 2, suffix='%(percent)d%%')
    bar.next()
    tcrrep.hcluster_df, tcrrep.Z = hcluster_diff(**default_hcluster_diff_kwargs)
    bar.next(); bar.finish()

    tcrrep.hcluster_df['prune'] = tcrrep.hcluster_df['K_neighbors'].apply(lambda x: 1 if (x < prune) else 0)#
    
    """ Do Basic Summary """
    mean_distance_      = list()
    percentage_node_25_ = list()
    percentage_node_50_ = list()
    percentage_node_75_ = list()
    n_rows = tcrrep.hcluster_df.shape[0]
    bar = IncrementalBar(f'Evaluate Clusters :', max = n_rows, suffix='%(percent)d%%')
    for i,r in tcrrep.hcluster_df.iterrows():
        bar.next()
        # <dfnode> is dataframe with all the clones at a given tree node
        dfnode = tcrrep.clone_df.iloc[r['neighbors_i'],]
        # <pwnod> is dataframe with all the clones at a given tree node
        pwnode = getattr(tcrrep, pwmat_str)[r['neighbors_i'],:][:,r['neighbors_i']]
        # get the non-diaganol entries.
        node_non_diag_entries = pwnode[~np.eye(pwnode.shape[0],dtype=bool)]
        # Compute the mean distance at the node 
        mean_distance_.append(str(round(node_non_diag_entries.mean(),1)))

        percentage_node_25 = 100*(node_non_diag_entries < 25).sum() / (node_non_diag_entries.size)
        percentage_node_50 = 100*(node_non_diag_entries < 50).sum() / (node_non_diag_entries.size)#100*((pwnode < 50).sum() - pwnode.shape[0]) / (pwnode.shape[0] * pwnode.shape[1])
        percentage_node_75 = 100*(node_non_diag_entries < 75).sum() / (node_non_diag_entries.size)#100*((pwnode < 100).sum() - pwnode.shape[0]) / (pwnode.shape[0] * pwnode.shape[1])
        percentage_node_25_.append(f"{round(percentage_node_25,1)}%")
        percentage_node_50_.append(f"{round(percentage_node_50,1)}%")
        percentage_node_75_.append(f"{round(percentage_node_75,1)}%")
    bar.next(); bar.finish()
    tcrrep.hcluster_df['mean_dist']   = mean_distance_
    tcrrep.hcluster_df['pct_dist_25'] =	percentage_node_25_
    tcrrep.hcluster_df['pct_dist_50'] =	percentage_node_50_
    tcrrep.hcluster_df['pct_dist_75'] =	percentage_node_75_
    
    """
    By default, treat each clone as a single entity if 'count_col' is single 
    """
    if default_member_summ_kwargs['count_col'] == 'single':
        single = True
        print("MAKING 'single' variable")
        tcrrep.hcluster_df['single'] = 1

    if single:
        tcrrep.clone_df['single'] = 1
        default_member_summ_kwargs['count_col'] = 'single'


    """ member_summ"""
    tcrrep.res_summary = \
        member_summ(res_df = tcrrep.hcluster_df,clone_df = tcrrep.clone_df, **default_member_summ_kwargs)


    tcrrep.hcluster_df_detailed = pd.concat([tcrrep.hcluster_df.copy(), tcrrep.res_summary.copy()], axis = 1)
    
    """ Add diversity stats"""
    tcrrep.clone_df['single'] = 1
    if single:
        fdiv75 = lambda ind : fuzzy_diversity(tcrrep.clone_df.iloc[ind,:]['single'], getattr(tcrrep, pwmat_str)[ind,:][:,ind],order=2, threshold=75)
        fdiv50  = lambda ind : fuzzy_diversity(tcrrep.clone_df.iloc[ind,:]['single'], getattr(tcrrep, pwmat_str)[ind,:][:,ind],order=2, threshold=50)
        fdiv25  = lambda ind : fuzzy_diversity(tcrrep.clone_df.iloc[ind,:]['single'], getattr(tcrrep, pwmat_str)[ind,:][:,ind],order=2, threshold=25)
    else:
        fdiv75 = lambda ind : fuzzy_diversity(tcrrep.clone_df.iloc[ind,:]['count'], getattr(tcrrep, pwmat_str)[ind,:][:,ind],order=2, threshold=75)
        fdiv50  = lambda ind : fuzzy_diversity(tcrrep.clone_df.iloc[ind,:]['count'], getattr(tcrrep, pwmat_str)[ind,:][:,ind],order=2, threshold=50)
        fdiv25  = lambda ind : fuzzy_diversity(tcrrep.clone_df.iloc[ind,:]['count'], getattr(tcrrep, pwmat_str)[ind,:][:,ind],order=2, threshold=25)
    tcrrep.hcluster_df_detailed['fuzzy_simpson_diversity_25']  = [str(round(fdiv25(ind),2)) for ind in tcrrep.hcluster_df_detailed.neighbors_i.to_list()]
    tcrrep.hcluster_df_detailed['fuzzy_simpson_diversity_50']  = [str(round(fdiv50(ind),2)) for ind in tcrrep.hcluster_df_detailed.neighbors_i.to_list()]
    tcrrep.hcluster_df_detailed['fuzzy_simpson_diversity_75'] = [str(round(fdiv75(ind),2)) for ind in tcrrep.hcluster_df_detailed.neighbors_i.to_list()]

    """Optional Add SVGs to hcluster_detailed"""
    if 'beta' in tcrrep.chains:
        _tcrsampler_svgs(tcrrep = tcrrep,
            default_background = None, 
            default_background_if_missing = None,
            cdr3_name = 'cdr3_b_aa',
            pwmat_str = pwmat_str_b, 
            chain = 'beta', 
            gene_names = ['v_b_gene','j_b_gene'],
            combine_olga = combine_olga)
    
    if 'alpha' in tcrrep.chains:
        _tcrsampler_svgs(tcrrep = tcrrep,
            default_background = None, 
            default_background_if_missing = None,
            cdr3_name = 'cdr3_a_aa',
            pwmat_str = pwmat_str_a, 
            chain = 'alpha', 
            gene_names = ['v_a_gene','j_a_gene'],
            combine_olga = combine_olga)


    """ Plot """
    html = plot_hclust_props(tcrrep.Z,
            res=tcrrep.hcluster_df_detailed,
            prune_col = 'prune',
            **default_plot_hclust_props)

    """ Write File """		
    with open(html_name, 'w') as fh:
        print(f"WRITING {html_name}")
        fh.write(html)


def report_kwargs(keyword_arguments):
    return ", ".join([f"{k} : {v}" for k,v in keyword_arguments.items()])

def _tcrsampler_svgs(
    tcrrep,
    default_background = None, 
    default_background_if_missing = None,
    cdr3_name = 'cdr3_b_aa',
    pwmat_str = 'pw_cdr3_b_aa', 
    chain = 'beta', 
    gene_names = ['v_b_gene','j_b_gene'],
    combine_olga = False,
    verbose = True):
    """
    Breath. What does this do?
    
    Given a TCRrep instance, this function samples a background repertoire
    using TCRsampler and makes svg-logos using palmotif. 

    This function doesn't return anything it. It needs to access 
    attribute values of a TCRrep (tcrrep) instance and 
    it modifies th etcrrep in place adding svgs and stats colums 
    to .hcluster_df_detailed DataFrame. TODO: could just output a dataframe 
    which would then just be concattenated.

    ONLY WORKS WITH _BETA using defaults:

    Notes
    -----
    Note: TCRSampler.build_background() accepts kwargs, we've set these as fixed as most user 
    won't know what these do and won't need to change them.
            max_rows : int
                Maximum clones per v,j pair (per subject)
            stratify_by_subject : bool
                If True, max_rows will apply to v,j,subject. If False, max_rows applies to v,j
            use_frequency : bool
                If True, uses frequency for ranking rows. If False, uses raw counts.
            make_singleton : bool
                If True, background is still sorted by frequency or counts, 
                but final fequency and counts values are overridden
                and set to 1.
    """  
    from tcrsampler.sampler import TCRsampler
    from palmotif import compute_pal_motif, svg_logo
    import pandas as pd
    from tcrdist.summarize import _select
    
    if chain == 'alpha' and tcrrep.organism == "mouse":
        # Here we enforce the rule that alpha-mouse cannot use an olga-sampler
        # TODO: This should be removed as soon as TCRsampler can be updated with a valid 
        # mouse-alpha simulated background.
        combine_olga = False
    
    # _default_sampler returns a TCRSampler based on organism and chain
    if verbose: print(f"INITIALIZING A TCRSAMPLER")
    print(tcrrep.organism, chain)
    t = _default_sampler(organism = tcrrep.organism, chain = chain)(
        default_background = default_background, 
        default_background_if_missing = default_background_if_missing)

    build_kargs = { 'max_rows' : 100,
                    'stratify_by_subject' : True,
                    'use_frequency' : True,
                    'make_singleton' : False}

    build_kargs_olga = {   'max_rows' : 1000,
                        'stratify_by_subject' : False,
                        'use_frequency' : True,
                        'make_singleton' : False}
    
    if verbose: print(f"BUILDING A DEEPER BACKGROUND {report_kwargs(build_kargs)}")

    t.build_background(**build_kargs)
    # Olga Sampler
    
    if combine_olga:
        t_olga = _default_sampler_olga(chain = chain, organism = tcrrep.organism)()
        if verbose: print(f"BUILDING A DEEPER BACKGROUND {report_kwargs(build_kargs_olga)}")
        t.build_background(**build_kargs_olga )

        olga_model = {
            ('beta','human')  : OlgaModel(recomb_type="VDJ", chain_folder = "human_T_beta"),
            ('alpha','human') : OlgaModel(recomb_type="VJ",  chain_folder = "human_T_alpha"),
            ('beta','mouse')  : OlgaModel(recomb_type="VDJ", chain_folder = "mouse_T_beta")
            }[(chain, tcrrep.organism)]



    if 'prune' not in tcrrep.hcluster_df.columns:
        if verbose: print("NO PRUNE COLUMNS USED ALL SET TO 0")
        tcrrep.hcluster_df['prune'] = 0
    
    print("ITERATE THROUGH CLUSTERS")
    svgs = list()
    svgs_raw = list()
    reference_unique = list()
    reference_unique_olga= list()
    reference_size = list()
    reference_size_olga = list()
    percent_missing_sampler = list()
    percent_missing_sampler_olga = list()
    n_rows = tcrrep.hcluster_df.shape[0]
    
    
    bar = IncrementalBar(f'Make {chain} SVGs :', max = n_rows, suffix='%(percent)d%%')
    for i,r in tcrrep.hcluster_df.iterrows():
        bar.next()
        if r['prune'] == 0:
            # <dfnode> is dataframe with all the clones at a given tree node
            dfnode   = tcrrep.clone_df.iloc[r['neighbors_i'],].copy()
            # <pwnode> Pairwise Matrix for node sequences
            pwnode   = getattr(tcrrep, pwmat_str)[r['neighbors_i'],:][:,r['neighbors_i']].copy()
            iloc_idx = pwnode.sum(axis = 0).argmin()
            centroid = dfnode[cdr3_name].to_list()[iloc_idx]
    
            # Compute gene usage at the node  
            # Convert to allele_01
            for gene_name in gene_names:
                dfnode[gene_name] = dfnode[gene_name].apply(lambda x : allele_01(x)) 
                      
            gene_usage = dfnode.groupby(gene_names).size() # e.g., ['v_b_gene','j_b_gene']
            gene_usage_tuples = gene_usage.reset_index().to_dict('split')['data']
            # Given gene usage use the <t> a TCRsampler instance to get background seqs
            
            # Adjust depth for small nodes
            adjust_depth = 10*round(10 / dfnode.shape[0])
            if adjust_depth < 10:
                adjust_depth = 10

            sampled_rep = t.sample( gene_usage_tuples,
                            flatten = True, depth = adjust_depth * 10)
            
            # Only keep the non-none sequences
            sampled_rep  = [x for x in sampled_rep if x is not None]
            # < missing_gene > Count the percentage missing, sampler returns none when no v,j pair is present
            expected_depth = dfnode.shape[0] * adjust_depth * 10
            recovered_depth = len(sampled_rep)
            percent_missing = round(100* (1- (recovered_depth / expected_depth)), 1)
                   
            percent_missing_sampler.append(f"{percent_missing}%")           
            reference_unique.append(str(pd.Series(sampled_rep).nunique()))
            reference_size.append(str(pd.Series(sampled_rep).count()))

            if combine_olga:
                # We modified Olga source code slightly, such that we simulated sequences 
                # with a given V,J gene usage
                # OLD METHOD WHERE WE ACTUALLY SAMPLED, slower but can go much deeper. I don't think one rare sequence however, really make a big difference.
                #flatten = lambda l: [item for sublist in l for item in sublist]
                #sampled_rep_olga = [olga_model.gen_cdr3s(allele_01(v),allele_01(j),n*adjust_depth*10) for v,j,n in gene_usage_tuples]
                #sampled_rep_olga = [x for x in flatten(sampled_rep_olga) if x is not None]
                sampled_rep_olga = t_olga.sample( gene_usage_tuples,
                            flatten = True, depth = adjust_depth * 10)

                sampled_rep_olga = [x for x in sampled_rep_olga if x is not None]
                
                expected_depth = dfnode.shape[0] * adjust_depth * 10
                recovered_depth = len(sampled_rep_olga)
                percent_missing_olga = round(100* (1- (recovered_depth / expected_depth)), 1)
                
                percent_missing_sampler_olga.append(f"{percent_missing_olga}%")            
                reference_unique_olga.append(str(pd.Series(sampled_rep_olga).nunique()))
                reference_size_olga.append(str(pd.Series(sampled_rep_olga).count()))
                
                # HERE WE COMBINE INTO A SINGLE BACKGROUND:
                sampled_rep = sampled_rep + sampled_rep_olga


            # Get motif matrix and motif stats
            motif, stat = compute_pal_motif(
                            seqs = _select(df = tcrrep.clone_df, 
                                        iloc_rows = r['neighbors_i'], 
                                        col = cdr3_name),
                            refs = sampled_rep, 
                            centroid = centroid)

            svgs.append(svg_logo(motif, return_str= True))

            # repeaat without references
            raw_motif, raw_stat = compute_pal_motif(
                seqs = _select(df = tcrrep.clone_df, 
                            iloc_rows = r['neighbors_i'], 
                            col = cdr3_name),
                centroid = centroid)
            # Convert the motif matrix into an svg_logo, append to list
            svgs_raw.append(svg_logo(raw_motif, return_str= True))
        else: 
            # If prune column is 1 don't go to the trouble of sampling and generating seqs
            svgs.append("PRUNE")
            svgs_raw.append("PRUNE")
            reference_size.append("PRUNE")
            reference_unique.append("PRUNE")
            percent_missing_sampler.append("PRUNE")
            percent_missing_sampler_olga.append("PRUNE")            
            reference_unique_olga.append("PRUNE") 
            reference_size_olga.append("PRUNE") 

    bar.next(); bar.finish()

    # The standard svg_ includes background, whereas raw has no background 
    tcrrep.hcluster_df_detailed[f'svg_{chain}'] 		= svgs
    tcrrep.hcluster_df_detailed[f'svg_raw_{chain}'] 	= svgs_raw
    tcrrep.hcluster_df_detailed[f'ref_size_{chain}'] 	= reference_size
    tcrrep.hcluster_df_detailed[f'ref_unique_{chain}'] 	= reference_unique
    tcrrep.hcluster_df_detailed[f'percent_missing_{chain}']  = percent_missing_sampler   
    if combine_olga:
        tcrrep.hcluster_df_detailed[f'ref_size_olga_{chain}'] 	      = reference_size_olga
        tcrrep.hcluster_df_detailed[f'ref_unique_olga_{chain}'] 	  = reference_unique_olga
        tcrrep.hcluster_df_detailed[f'percent_missing_olga_{chain}']  = percent_missing_sampler_olga  
    
    return True



def _default_sampler(organism = 'human', chain = 'beta'):
    assert organism in ['human', 'mouse']
    assert chain in ['beta','alpha']

    default_tcrsampler_generator = {
        ('human','beta'): 
            _default_tcrsampler_human_beta,
        ('human','alpha'): 
            _default_tcrsampler_human_alpha,
        ('mouse','beta'): 
            _default_tcrsampler_mouse_beta,
        ('mouse','alpha'): 
            _default_tcrsampler_mouse_alpha, 
        }[(organism, chain)]
    
    return default_tcrsampler_generator

def _default_sampler_olga(organism = 'human', chain = 'beta'):
    assert organism in ['human', 'mouse']
    assert chain in ['beta','alpha']
    if organism == "mouse" and chain == "alpha":
        raise KeyError("No currenlty available mouse-alpha chain background form olga")

    default_tcrsampler_generator = {
        ('human','beta'): 
            _default_tcrsampler_olga_human_beta,
        ('human','alpha'): 
            _default_tcrsampler_olga_human_alpha,
        ('mouse','beta'): 
            _default_tcrsampler_olga_mouse_beta
        }[(organism, chain)]
    
    return default_tcrsampler_generator


def _default_tcrsampler_olga_human_beta(default_background = None, default_background_if_missing=None):
    """
    Responsible for providing the default human beta sampler 'britanova_human_beta_t_cb.tsv.sampler.tsv'

    Returns
    -------
    t : tcrsampler.sampler.TCRsampler 
    """
    from tcrsampler.sampler import TCRsampler
    if default_background is None:
        default_background =  'olga_human_beta_t.sampler.tsv'
        
    if default_background_if_missing is None:
        default_background_if_missing ='olga_sampler.zip'
    
    print(default_background)

    try: 
        t = TCRsampler(default_background=default_background)
    except OSError:
        t = TCRsampler()
        t.download_background_file(default_background_if_missing)
        t = TCRsampler(default_background=default_background)
    return t

def _default_tcrsampler_olga_human_alpha(default_background = None, default_background_if_missing=None):
    """
    Responsible for providing the default human beta sampler 'britanova_human_beta_t_cb.tsv.sampler.tsv'

    Returns
    -------
    t : tcrsampler.sampler.TCRsampler 
    """
    from tcrsampler.sampler import TCRsampler
    if default_background is None:
        default_background =  'olga_human_alpha_t.sampler.tsv'
        
    if default_background_if_missing is None:
        default_background_if_missing ='olga_sampler.zip'
    
    
    print(default_background)

    try: 
        t = TCRsampler(default_background=default_background)
    except OSError:
        t = TCRsampler()
        t.download_background_file(default_background_if_missing)
        t = TCRsampler(default_background=default_background)
    return t

def _default_tcrsampler_olga_mouse_beta(default_background = None, default_background_if_missing=None):
    """
    Responsible for providing the default human beta sampler 'britanova_human_beta_t_cb.tsv.sampler.tsv'

    Returns
    -------
    t : tcrsampler.sampler.TCRsampler 
    """
    from tcrsampler.sampler import TCRsampler
    if default_background is None:
        default_background =  'olga_mouse_beta_t.sampler.tsv'
        
    if default_background_if_missing is None:
        default_background_if_missing ='olga_sampler.zip'
    
    
    print(default_background)

    try: 
        t = TCRsampler(default_background=default_background)
    except OSError:
        t = TCRsampler()
        t.download_background_file(default_background_if_missing)
        t = TCRsampler(default_background=default_background)
    return t


def _default_tcrsampler_human_beta(default_background = None, default_background_if_missing=None):
    """
    Responsible for providing the default human beta sampler 'britanova_human_beta_t_cb.tsv.sampler.tsv'

    Returns
    -------
    t : tcrsampler.sampler.TCRsampler 
    """
    from tcrsampler.sampler import TCRsampler
    if default_background is None:
        default_background = 'britanova_human_beta_t_cb.tsv.sampler.tsv'
        
    if default_background_if_missing is None:
        default_background_if_missing ='britanova_human_beta_t_cb.tsv.sampler.tsv.zip'
    
    
    print(default_background)

    try: 
        t = TCRsampler(default_background=default_background)
    except OSError:
        t = TCRsampler()
        t.download_background_file(default_background_if_missing)
        t = TCRsampler(default_background=default_background)
    return t

def _default_tcrsampler_human_alpha(default_background = None, default_background_if_missing=None ):
    """
    Responsible for providing the default human alpha sampler 'ruggiero_human_alpha_t.tsv.sampler.tsv'
    """
    from tcrsampler.sampler import TCRsampler
    if default_background is None:
        default_background = 'ruggiero_human_alpha_t.tsv.sampler.tsv'
    if default_background_if_missing is None:
        default_background_if_missing =  'ruggiero_human_alpha_t.tsv.sampler.tsv.zip'
    
    print(default_background)

    try: 
        t = TCRsampler(default_background=default_background)
    except OSError:
        t = TCRsampler()
        t.download_background_file(default_background_if_missing)
        t = TCRsampler(default_background=default_background)
    return t

def _default_tcrsampler_mouse_beta(default_background = None, default_background_if_missing=None):
    """
    Responsible for providing the default mouse beta sampler

    Returns
    -------
    t : tcrsampler.sampler.TCRsampler 
    """
    from tcrsampler.sampler import TCRsampler

    if default_background is None:
        default_background = 'ruggiero_mouse_beta_t.tsv.sampler.tsv'
    if default_background_if_missing is None:
        default_background_if_missing =  'ruggiero_mouse_sampler.zip'

    print(default_background)

    try: 
        t = TCRsampler(default_background=default_background)
    except OSError:
        t = TCRsampler()
        t.download_background_file(default_background_if_missing)
        t = TCRsampler(default_background=default_background)
    return t

def _default_tcrsampler_mouse_alpha(default_background = None, default_background_if_missing=None):
    """
    Responsible for providing the default mouse alpha sampler
    """
    from tcrsampler.sampler import TCRsampler
    
    if default_background is None:
        default_background = 'ruggiero_mouse_alpha_t.tsv.sampler.tsv'
    if default_background_if_missing is None:
        default_background_if_missing =  'ruggiero_mouse_sampler.zip'
    
    print(default_background)

    try: 
        t = TCRsampler(default_background=default_background)
    except OSError:
        t = TCRsampler()
        t.download_background_file(default_background_if_missing)
        t = TCRsampler(default_background=default_background)
    return t

def allele_01(genename):
    """
    >>> allele_01('TRBV19*02')
    'TRBV19*01'
    """
    g,a = re.match(pattern = '(.*)([0-9])$', string= genename).groups()
    return(f"{g}1")