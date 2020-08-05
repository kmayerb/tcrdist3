import pandas as pd
import pwseqdist as pw
import numpy as np
import warnings
from . import repertoire_db
from tcrdist.rep_funcs import _pws
from zipdist.zip2 import Zipdist2
import sys

class TCRrep:
    """
    Flexible distance measures for comparing T cell receptors

    The TCRrep Class hold T cell repertoire data, infers CDRs from v-gene name, and 
    computes multi-CDR 'tcrdistance'.  

    Attributes
    ----------
    cell_df : pd.DataFrame or None
        Pandas DataFrame containing cell level information 
    clone_df : pd.DataFrame or None
        Pandas DataFrame containing clone level information. 
        This can be provided directly from a program like MIXCR or
        can be inferred by deduplicating a cell_df. 
    organism = str, 
        specifies relevant organism for analysis: 'human' or 'mouse'
    chains : list
        specifies relevant chains for single or paried TCR analysis
        ['alpha','beta'], ['alpha'], ['beta'], ['gamma','delta'],  ['gamma'] or ['delta']  
    db_file : str
        specifies refereence file. The default is 'alphabeta_gammadelta_db.tsv' which
        is preinstalled  with the install python3.7/site-packages/tcrdist/db/
    archive_name : str
        Name for archive file. (only used if archive result is True)
    archive_result: bool
        If True save result to .tar.gz archive
    imgt_aligned : bool
        If True, by default, cdr1, cdr2,and pmhc are inferred aligned to fixed length with gaps.
        If False, cdr1, cdr2,and pmhc are returned as ungapped raw sequences.   
    infer_all_genes : bool
        If True, load all_gene reference from 'db_file`.
    infer_cdrs : bool
        If True, infer cdr1, cdr2,and pmhc from the v gene name
    infer_index_cols : bool
        If True, infer index_cols used to deduplicate cell_df.
        If False, index_cols can be specified directly after initialization.
    index_cols : list
        list of index colums used to deduplicate cell_df to clone_df
    deduplicate : bool
        If True, then clone_df will be assigned cell_df grouped_by
        index_cols.
        If False, and clone_df is None, then clone_df will be be 
        assigned a copy of cell_df.
    use_defaults : True
        If True, use default metrics, weights, and kargs
    store_all_cdr : True,
        If True, save all distance matrices for each CDR (e.g., pw_cdr3_b_aa).
        If False, only save pw_alpha and pw_beta
    compute_distances : True
        If True, automatically compute distances
    cpus : int, 
        Number of cpus to use. In general 1 cpu is sufficient from default Numba metric 
        with less than 10^7  pairwise compairsons. However, more cpus will 
        result in a speedup for metrics like pw.metrics.nw_hamming_metric for more than 10^6 
        pairwise compairsons.
    
    Example
    -------
    .. code-block:: python

        import pandas as pd
        from tcrdist.repertoire import TCRrep

        df = pd.read_csv("dash.csv")
        tr = TCRrep(
                cell_df = df, 
                organism = 'mouse', 
                chains = ['alpha','beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv',
                clone_df          = None,
                imgt_aligned      = True,
                infer_all_genes   = True,
                infer_cdrs        = True,
                infer_index_cols  = True,
                deduplicate       = True,
                use_defaults      = True,
                store_all_cdr     = True,
                compute_distances = True,
                index_cols        = None,
                cpus              = 1, 
                db_file           = 'alphabeta_gammadelta_db.tsv')

    Notes
    -----
    The default initialization compute pairwise distance matrices. 

    See examples at https://tcrdist3.readthedocs.io/ for more information
    about flexibility of TCRrep.

    Computed pairwised distance matrics that can be accessed:
        TCRrep.pw_alpha, 
        TCRrep.pw_beta,  
        TCRrep.pw_cdr3_a_aa, 
        TCRrep.pw_cdr3_b_aa  
    """
    def __init__(self,
                 organism          = "mouse",
                 chains            = ['alpha', 'beta'],
                 db_file           = 'alphabeta_gammadelta_db.tsv',
                 archive_name      = 'tcrdist3.archive',
                 blank             = False,
                 cell_df           = None,
                 clone_df          = None,
                 imgt_aligned      = True,
                 infer_all_genes   = True,
                 infer_cdrs        = True,
                 infer_index_cols  = True,
                 deduplicate       = True,
                 use_defaults      = True,
                 store_all_cdr     = True,
                 compute_distances = True,
                 index_cols        = None,
                 cpus              = 1,
                 df2               = None, 
                 archive_result    = False):
        

        self.organism = organism
        self._validate_organism()   

        self.chains = chains
        self._validate_chains()
        
        self.archive_name = archive_name
        self.blank = blank 
        # blank is only used if reloading from .tar.gz
        if self.blank:
            self.cell_df           = None
            self.clone_df          = None
            self.imgt_aligned      = False
            self.infer_all_genes   = False
            self.store_all_cdr     = False
            self.index_cols        = None
            self.cpus              = 1,
            self.df2               = None 
            self.archive_result    = False
        else:
            self.index_cols = index_cols
            self.clone_df = clone_df 

            if cell_df is None:
                cell_df = pd.DataFrame()
            self.cell_df = cell_df
            self._validate_cell_df()
            self.clone_df = clone_df
            self.df2 = df2
            
            
            self.db_file = db_file
            self._validate_db_file()

            self.imgt_aligned = imgt_aligned
            self._validate_imgt_aligned()

            self.store_all_cdr = store_all_cdr
            self.cpus = cpus

            self.archive_result = archive_result
            self.infer_all_genes = infer_all_genes
            self.use_defaults = use_defaults

            if infer_all_genes:     
                self.all_genes = repertoire_db.RefGeneSet(db_file).all_genes
            
            if infer_cdrs:
                for chain in self.chains:
                    self.infer_cdrs_from_v_gene(chain = chain, imgt_aligned = self.imgt_aligned)
                        # Assume all provided columns are index columns, except 'count' 'cell_id', 'clone_id'
            
            if infer_index_cols:
                self.infer_index_cols()
            
            if deduplicate:
                self.infer_index_cols()
                self.deduplicate()
            else: 
                if self.clone_df is None:
                    self.clone_df = self.cell_df.copy()

            if use_defaults:
                self._initialize_chain_specific_attributes()
        
            if compute_distances:
                self.compute_distances()

            if self.archive_result:
                self.archive()
     
    def archive(self, 
                dest = None, 
                dest_tar_name = None, 
                verbose = True, 
                use_csv = False):
        """ 
        Use Zipdist2 to Make an Archive.tar.gz 

        Parameters
        ----------
        dest : str
            e.g., 'default_archive'
        dest_tar_name : str
            e.g., 'default_archive.tar.gz'
        verbose : bool
            if True, report steps in archive process
        use_csv : bool
            if True, archive will include .csv file. Useful for porting files to other applications, but creates large files.

        Example
        -------
        .. code-block:: python
            tr = TCRrep(cell_df = pd.DataFrame(), organism = "mouse")
            tr.archive(dest = "default_archive", dest_tar_name = "default_archive.tar.gz")

        Notes
        -----
        See :py:meth:`tcrdist.repertoire.rebuild`: for reubilding a TCRrep instance from 
        an TCRrep archive .tar.gz file.
        
        """
        if dest is None:
            dest = self.archive_name
        if dest_tar_name is None:
            dest_tar_name = f"{dest}.tar.gz"

        self.cell_df_index = self.cell_df.index.copy()
        self.cell_df = self.cell_df.reset_index()

        z = Zipdist2(name = dest_tar_name , target = self)
        z._save(dest = dest, dest_tar = dest_tar_name, verbose = verbose, use_csv = use_csv )
        sys.stdout.write(f"\tArchiving your TCRrep using Zipdist2 in [{dest_tar_name}]\n")
    
    def rebuild(self,
                dest = None,
                dest_tar_name = None, 
                verbose = True ):
        """ 
        Use Zipdist2 to reubild a TCRrep instance from an Archive.tar.gz

        Parameters
        ----------
        dest_tar_name : str
            e.g., 'default_archive.tar.gz'
        verbose : bool
            If True, report rebuilding process steps.

        Example
        -------
        Shows :py:meth:`tcrdist.repertoire.archive` and :py:meth:`tcrdist.repertoire.rebuild` 
        used together.

        .. code-block:: python
        
            tr = TCRrep(cell_df = pd.DataFrame(), organism = "mouse")
            tr.archive(dest = "default_archive", dest_tar_name = "default_archive.tar.gz")
            tr_new = TCRrep(cell_df = pd.DataFrame(), organism = "mouse")
            tr_new.rebuild(dest_tar_name = "default_archive.tar.gz")
        
        Notes
        -----
        See :py:meth:`tcrdist.repertoire.archive` for creating TCRrep archive file.
        """
        #tr = TCRrep(cell_df=df.iloc[0:0,:], chains=chains, organism='mouse')
        if dest is None:
            dest = self.archive_name
        if dest_tar_name is None:
            dest_tar_name = f"{dest}.tar.gz"

        z = Zipdist2(name = dest, target = self)
        z._build(dest_tar = dest_tar_name , target = self, verbose = verbose)
        
        # VALIDATION OF INPUTS
        # check that chains are valid.
        self._validate_organism()
        self._validate_chains()
        # check that  is a pd.DataFrame
        self._validate_cell_df()

        # RE INIT the REFERENCE DB see repertoire_db.py
        self.all_genes = repertoire_db.RefGeneSet(self.db_file).all_genes
       

    def compute_distances(self, df = None):
        if df is None:
            df = self.clone_df
        if 'alpha' in self.chains:
            pw_alpha  = _pws(
                df = df,
                metrics = self.metrics_a, 
                weights = self.weights_a, 
                kargs   = self.kargs_a, 
                cpu     = self.cpus, 
                store   = self.store_all_cdr)
            self._assign_distance_attributes(d = pw_alpha, chain = 'alpha')
        if 'beta' in self.chains:
            pw_beta = _pws(
                df = df, 
                metrics = self.metrics_b, 
                weights = self.weights_b, 
                kargs   = self.kargs_b, 
                cpu     = self.cpus, 
                store   = self.store_all_cdr)
            self._assign_distance_attributes(d = pw_beta, chain = 'beta')
        if 'gamma' in self.chains:
            pw_gamma = _pws(
                df = df,
                metrics = self.metrics_g, 
                weights = self.weights_g, 
                kargs   = self.kargs_g, 
                cpu     = self.cpus, 
                store   = self.store_all_cdr)
            self._assign_distance_attributes(d = pw_gamma, chain = 'gamma')
        if 'delta' in self.chains:
            pw_delta = _pws(
                df = df, 
                metrics = self.metrics_d,
                weights = self.weights_d, 
                kargs   = self.kargs_d, 
                cpu     = self.cpus, 
                store   = self.store_all_cdr)
            self._assign_distance_attributes(d = pw_delta, chain = 'delta')
    
    def compute_rect_distances(self, df = None, df2 = None, store = None):
        if df is None:
            df = self.clone_df
        if df2 is None:
            df2 = self.df2
        if store is None:
            store = self.store_all_cdr
        if 'alpha' in self.chains:
            pw_alpha  = _pws(
                df = df,
                df2 = df2,
                metrics = self.metrics_a, 
                weights = self.weights_a, 
                kargs   = self.kargs_a, 
                cpu     = self.cpus, 
                store   = store)
            self._assign_distance_attributes(d = pw_alpha, chain = 'alpha',prefix = "rw")
        if 'beta' in self.chains:
            pw_beta = _pws(
                df = df,
                df2 = df2, 
                metrics = self.metrics_b, 
                weights = self.weights_b, 
                kargs   = self.kargs_b, 
                cpu     = self.cpus, 
                store   = store)
            self._assign_distance_attributes(d = pw_beta, chain = 'beta',prefix = "rw")
        if 'gamma' in self.chains:
            pw_gamma = _pws(
                df = df,
                df2 = df2, 
                metrics = self.metrics_g, 
                weights = self.weights_g, 
                kargs   = self.kargs_g, 
                cpu     = self.cpus, 
                store   = store)
            self._assign_distance_attributes(d = pw_gamma, chain = 'gamma', prefix = "rw")
        if 'delta' in self.chains:
            pw_delta = _pws(
                df = df,
                df2 = df2,
                metrics = self.metrics_d,
                weights = self.weights_d, 
                kargs   = self.kargs_d, 
                cpu     = self.cpus, 
                store   = store)
            self._assign_distance_attributes(d = pw_delta, chain = 'delta', prefix = "rw")

    def _assign_distance_attributes(self, d:dict, chain:str, prefix = 'pw'):
        for k in d.keys():
            if k == 'tcrdist':
                pw_k = f"{prefix}_{chain}"
            else:
                pw_k = f"{prefix}_{k}"
            setattr(self, pw_k, d[k])

    def infer_cdrs_from_v_gene(self, chain, imgt_aligned = True):
        """
        Function taking TCR v-gene name to infer the amino amino_acid
        sequence of cdr1, cdr2, and pmhc loop regions.

        Parameters
    	----------
        chain : string
            'alpha', 'beta', 'gamma', or 'delta'
        imgt_aligned : boolean
            if True cdr1, cdr2, cdr2.5 will be returned with gaps
            and by definition will be the same length. MSH.......ET

        Returns
    	-------
        self.cell_df : pandas.core.frame.DataFrame
    	   Assigns [cdr3|cdr2|cdr1|pmhc]_[a|b|d|g]_aa columns in TCRrep.cell_df

        Examples
    	--------
        >>> testrep = TCRrep(cell_df = example_df, organism = "human", chains= ["alpha","beta"])
        >>> testrep.infer_cdrs_from_v_gene(chain = "alpha")
        >>> testrep.infer_cdrs_from_v_gene(chain = "beta")

        Notes
    	-----
        This function takes the v-gene names and infers the amino acid
        sequence of the cdr1, cdr2, and pmhc region (pmhc refers to the
        pMHC-facing loop between CDR2 and CDR3 (IMGT alignment columns 81 - 86.
      
        These sequences are based up on lookup from the dictionary here:

        self.generate_ref_genes_from_db(db_file)

        imgt_aligned : boolean
            if True cdr1, cdr2, cdr2.5 will be returned with gaps
            and by definition will be the same length.
            MSH.......ET
            FNH.......DT
            LGH.......NA

        References
        ----------

        IMGT definitions of cdr1, cdr2, and pMHC-facing can be found here
        http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
        """

        if not imgt_aligned:
            self.imgt_aligned_status = False
            f0 = lambda v : self._map_gene_to_reference_seq2(gene = v,
                                                             cdr = 0,
                                                             organism = self.organism,
                                                             attr ='cdrs_no_gaps')
            f1 = lambda v : self._map_gene_to_reference_seq2(gene = v,
                                                             cdr = 1,
                                                             organism = self.organism,
                                                             attr ='cdrs_no_gaps')
            f2 = lambda v : self._map_gene_to_reference_seq2(gene = v,
                                                             cdr = 2,
                                                             organism = self.organism,
                                                             attr ='cdrs_no_gaps')
        else:
            self.imgt_aligned_status = True
            f0 = lambda v : self._map_gene_to_reference_seq2(gene = v,
                                                             cdr = 0,
                                                             organism = self.organism,
                                                             attr ='cdrs')
            f1 = lambda v : self._map_gene_to_reference_seq2(gene = v,
                                                             cdr = 1,
                                                             organism = self.organism,
                                                             attr ='cdrs')
            f2 = lambda v : self._map_gene_to_reference_seq2(gene = v,
                                                             cdr = 2,
                                                             organism = self.organism,
                                                             attr ='cdrs')
        if chain == "alpha":
            self.cell_df['cdr1_a_aa'] = list(map(f0, self.cell_df.v_a_gene))
            self.cell_df['cdr2_a_aa'] = list(map(f1, self.cell_df.v_a_gene))
            self.cell_df['pmhc_a_aa'] = list(map(f2, self.cell_df.v_a_gene))
        if chain == "beta":
            self.cell_df['cdr1_b_aa'] = list(map(f0, self.cell_df.v_b_gene))
            self.cell_df['cdr2_b_aa'] = list(map(f1, self.cell_df.v_b_gene))
            self.cell_df['pmhc_b_aa'] = list(map(f2, self.cell_df.v_b_gene))
        if chain == "gamma":
            self.cell_df['cdr1_g_aa'] = list(map(f0, self.cell_df.v_g_gene))
            self.cell_df['cdr2_g_aa'] = list(map(f1, self.cell_df.v_g_gene))
            self.cell_df['pmhc_g_aa'] = list(map(f2, self.cell_df.v_g_gene))
        if chain == "delta":
            self.cell_df['cdr1_d_aa'] = list(map(f0, self.cell_df.v_d_gene))
            self.cell_df['cdr2_d_aa'] = list(map(f1, self.cell_df.v_d_gene))
            self.cell_df['pmhc_d_aa'] = list(map(f2, self.cell_df.v_d_gene))
    
    def infer_index_cols(self):
        """
        Infers index columns from TCRrep.cell_df
        """
        self.index_cols = [item for item in self.cell_df.columns.to_list() if item not in ['count', 'cell_id', 'clone_id']]
        
    def deduplicate(self):
        """
        With attribute self.index_col calls _deduplicate() and assigns result to attribute self.clone_df
        """
        clone_df = _deduplicate(self.cell_df, self.index_cols)
        
        # check if any clones were lost due to missing information
        if np.sum(self.cell_df['count']) != np.sum(clone_df['count']):
            n_cells_lost = np.sum(self.cell_df['count']) - np.sum(clone_df['count'])
            n_cell = np.sum(self.cell_df['count'])
            warnings.warn(f"Not all cells/sequences could be grouped into clones. {n_cells_lost} of {n_cell} were not captured. This occurs when any of the values in the index columns are null or missing for a given sequence. To see entries with missing values use: tcrdist.repertoire.TCRrep._show_incomplete()\n")
        
        # if no clone id column provided thetrn create one as a sequence of numbers
        if "clone_id" not in clone_df:
            N = clone_df.shape[0]
            clone_df['clone_id'] = range(1, N + 1 ,1)
        
        self.clone_df = clone_df
        return clone_df.copy()

    def show_incomplete(self):
        """
        Returns a dataframe with those cells that do not have a valid entry 
        for any one of the specified index columns
        """   
        ind = self.cell_df[self.index_cols].isnull().any(axis = 1)   
        incomplete_clones = self.cell_df.loc[ind,self.index_cols].copy()
        return incomplete_clones 

    """
        _INTERNAL functions
    """

    def _initialize_chain_specific_attributes(self):
        """
        Initialize pw object and default substitution matrix (smat) based on
        chains arguments.

        Naming of all objects have a standardized order
            region_chain_molecular_object
            (cdr3)_(a|b|d|g)_(aa|p)_(pw|smat|hmat)

        """
        if "alpha" in self.chains:
            self.cdr3_a_aa_smat = 'blosum62'
            self.cdr2_a_aa_smat = 'blosum62'
            self.cdr1_a_aa_smat = 'blosum62'
            self.pmhc_a_aa_smat = 'blosum62'
            self.index_cols.append("cdr3_a_aa")
            self.metrics_a = { 
                "cdr3_a_aa" : pw.metrics.nb_vector_tcrdist,
                "pmhc_a_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr2_a_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr1_a_aa" : pw.metrics.nb_vector_tcrdist }
            self.weights_a = { 
                "cdr3_a_aa" : 3,
                "pmhc_a_aa" : 1,
                "cdr2_a_aa" : 1,
                "cdr1_a_aa" : 1}
            self.kargs_a = {
                "cdr3_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':3, 'ctrim':2, 'fixed_gappos':False},
                "pmhc_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr2_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr1_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True}}
        
        if 'beta' in self.chains:
            self.cdr3_b_aa_smat = 'blosum62'
            self.cdr2_b_aa_smat = 'blosum62'
            self.cdr1_b_aa_smat = 'blosum62'
            self.pmhc_b_aa_smat = 'blosum62'
            self.index_cols.append("cdr3_b_aa")
            self.metrics_b = { "cdr3_b_aa" : pw.metrics.nb_vector_tcrdist,
                               "pmhc_b_aa" : pw.metrics.nb_vector_tcrdist,
                               "cdr2_b_aa" : pw.metrics.nb_vector_tcrdist,
                               "cdr1_b_aa" : pw.metrics.nb_vector_tcrdist }
            self.metrics_b = { 
                "cdr3_b_aa" : pw.metrics.nb_vector_tcrdist,
                "pmhc_b_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr2_b_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr1_b_aa" : pw.metrics.nb_vector_tcrdist }
            self.weights_b = { 
                "cdr3_b_aa" : 3,
                "pmhc_b_aa" : 1,
                "cdr2_b_aa" : 1,
                "cdr1_b_aa" : 1}
            self.kargs_b = {
                "cdr3_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':3, 'ctrim':2, 'fixed_gappos':False},
                "pmhc_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr2_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr1_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True}}

        if 'gamma' in self.chains:
            self.cdr3_g_aa_smat = 'blosum62'
            self.cdr2_g_aa_smat = 'blosum62'
            self.cdr1_g_aa_smat = 'blosum62'
            self.pmhc_g_aa_smat = 'blosum62'
            self.index_cols.append("cdr3_g_aa")
            self.metrics_g = { "cdr3_g_aa" : pw.metrics.nb_vector_tcrdist,
                               "pmhc_g_aa" : pw.metrics.nb_vector_tcrdist,
                               "cdr2_g_aa" : pw.metrics.nb_vector_tcrdist,
                               "cdr1_g_aa" : pw.metrics.nb_vector_tcrdist }
            self.metrics_g = { 
                "cdr3_g_aa" : pw.metrics.nb_vector_tcrdist,
                "pmhc_g_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr2_g_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr1_g_aa" : pw.metrics.nb_vector_tcrdist }
            self.weights_g = { 
                "cdr3_g_aa" : 3,
                "pmhc_g_aa" : 1,
                "cdr2_g_aa" : 1,
                "cdr1_g_aa" : 1}
            self.kargs_g = {
                "cdr3_g_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':3, 'ctrim':2, 'fixed_gappos':False},
                "pmhc_g_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr2_g_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr1_g_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True}}
      

        if 'delta' in self.chains:
            self.cdr3_d_aa_smat = 'blosum62'
            self.cdr2_d_aa_smat = 'blosum62'
            self.cdr1_d_aa_smat = 'blosum62'
            self.pmhc_d_aa_smat = 'blosum62'
            self.index_cols.append("cdr3_d_aa")
            self.metrics_d = { "cdr3_d_aa" : pw.metrics.nb_vector_tcrdist,
                               "pmhc_d_aa" : pw.metrics.nb_vector_tcrdist,
                               "cdr2_d_aa" : pw.metrics.nb_vector_tcrdist,
                               "cdr1_d_aa" : pw.metrics.nb_vector_tcrdist }
            self.metrics_d = { 
                "cdr3_d_aa" : pw.metrics.nb_vector_tcrdist,
                "pmhc_d_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr2_d_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr1_d_aa" : pw.metrics.nb_vector_tcrdist }
            self.weights_d = { 
                "cdr3_d_aa" : 3,
                "pmhc_d_aa" : 1,
                "cdr2_d_aa" : 1,
                "cdr1_d_aa" : 1}
            self.kargs_d = {
                "cdr3_d_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':3, 'ctrim':2, 'fixed_gappos':False},
                "pmhc_d_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr2_d_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr1_d_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True}}
      

    def _map_gene_to_reference_seq2(self,
                                    organism,
                                    gene,
                                    cdr,
                                    attr = 'cdrs_no_gaps'):
        """
        Internal function that looks up the cdr sequence (gapped or ungapped)
        from the self.all_genes library

        Parameter
        ---------
        organism : string
            mouse or human
        gene : string
            specifies the TCR gene such as 'TRAV1*01'
        cdr : int
            0 - CDR1, 1-CDR2 and 2 - CDR2.5
        attr : string
            'cdrs_no_gaps' or 'cdrs_aligned' with gaps from IMGT
        """
        try:
            aa_string = self.all_genes[organism][gene].__dict__[attr][cdr]
        except KeyError:
            aa_string = None
            warnings.warn("{} gene was not recognized in reference db no cdr seq could be inferred".format(gene))
        return(aa_string)
       
        """
        _VALIDATION functions - check validity of inputs and communicate errors to user
        """

    def _validate_db_file(self):
        """
        Issues warning if invalid organism is passed to TCRrep __init__
        """
        if self.db_file not in ['alphabeta_gammadelta_db.tsv','alphabeta_db.tsv','gammadelta_db.tsv']:
            warnings.warn("db_file must be 'alphabeta_gammadelta_db.tsv' or 'alphabeta_db.tsv' or 'gammadelta_db.tsv' unless you have built tcrdist3 from scratch")
    
    def _validate_organism(self):
        """
        Raises 
        ------
        ValueError if invalid organism is passed to TCRrep __init__
        """
        if self.organism not in ["human", "mouse"]:
            raise ValueError("organism must be 'mouse' or 'human'")
    
    def _validate_chains(self):
        """
        Raises
        ------ 
        ValueError if invalid chains are passed to TCRrep __init__
        """
        check_chains_arg = ['alpha', 'beta', "gamma", "delta"]
        if len([c for c in self.chains if c not in check_chains_arg]) > 0:
            raise ValueError('TCRrep chains arg can be one or more of the '
                                'following {} case-sensitive'.format(check_chains_arg))
    def _validate_cell_df(self):
        """
        Issues warning if TCRrep.cell_df not properly formatted.
        """
        if not isinstance(self.cell_df, pd.DataFrame):
            warnings.warn('TCRrep cell_df argument must be pandas.DataFrame unless you are providing a clone_df directly\n')

        else:
            cell_df_columns = self.cell_df.columns.to_list()
            
            if "count" not in cell_df_columns:
                warnings.warn("cell_df needs a counts column to track clonal number of frequency\n")
            
            if "alpha" in self.chains:
                if not "cdr3_a_aa" in cell_df_columns:
                    warnings.warn("cell_df needs a column called 'cdr3_a_aa' to track the CDR3 amino acid sequence\n")
                if not "v_a_gene" in cell_df_columns:
                    warnings.warn("cell_df needs a column called 'v_a_gene' for default functions\n")
            if "beta" in self.chains:
                if not "cdr3_b_aa" in cell_df_columns:
                    warnings.warn("cell_df needs a column called 'cdr3_b_aa' to track the CDR3 amino acid sequence\n")
                if not "v_b_gene" in cell_df_columns:
                    warnings.warn("cell_df needs a column called 'v_b_gene' for default functions\n")
            if "gamma" in self.chains:
                if not "cdr3_g_aa" in cell_df_columns:
                    warnings.warn("cell_df needs a column called 'cdr3_g_aa' to track the CDR3 amino acid sequence\n")
                if not "v_g_gene" in cell_df_columns:
                    warnings.warn("cell_df needs a column called 'v_g_gene' for default functions\n")
            if "delta" in self.chains:
                if not "cdr3_d_aa" in cell_df_columns:
                    warnings.warn("cell_df needs a column called 'cdr3_d_aa' to track the CDR3 amino acid sequence\n")
                if not "v_d_gene" in cell_df_columns:
                    warnings.warn("cell_df needs a column called 'v_d_gene' for default functions\n")

    def _validate_imgt_aligned(self):
        """
        Raises
        ------
        ValueError if TCRrep.imgt_aligned is not a boolean
        """
        if not isinstance(self.imgt_aligned, bool):
            raise ValueError('TCRrep imgt_aligned argument must be a boolean')

    # def _validate_cdr_sequences():
    #     pass
    #     # TODO: Avoid trouble by wanring the user that they have non valid characters 

    
def _deduplicate(cell_df, index_cols):
    """
    Use index_cols to group by and group identical entries. The input DataFrame
    must have a column 'count'.
    """
    clones = cell_df.groupby(index_cols)['count'].agg(np.sum).reset_index()
    return clones