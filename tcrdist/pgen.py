import pandas as pd
import numpy as np
import os.path as op
import os
import olga
import olga.load_model as load_model
import olga.generation_probability as generation_probability
import olga.generation_probability as pgen
from olga.utils import nt2aa, determine_seq_type
from tcrdist.paths import path_to_default_models
path_to_olga_default_models = path_to_default_models

class OlgaModel:
    """
    A Class One Layer Up From Olga Classes To Help Store Configuration

    Attributes
    ----------
    chain_folder : string
        name of chain specific folder in the olga_default_models folder
    recomb_type : string
            'VDJ' or "VJ"
    """


    def __init__(self, chain_folder, recomb_type):
        """
        Sets up an OlgaModel that can be used multiple times. For instance to
        generate generation probabilities for 10K sequences.

        chain_folder : string
            'human_T_beta', 'human_T_alpha'
        recomb_type : string
            'VDJ' or "VJ"
        """

        self.chain_folder     = chain_folder
        self.recomb_type      = recomb_type
        self.generative_model = None
        self.genomic_data     = None
        self.pgen_model       = None


        self._validate_chain_folder_arg()
        self._validate_recomb_type_arg()
        self._validate_chain_folder_with_recomb_type()



        params_file_name = op.join(path_to_olga_default_models,
                                   chain_folder,
                                   'model_params.txt')
        marginals_file_name = op.join(path_to_olga_default_models,
                                      chain_folder,
                                      'model_marginals.txt')
        V_anchor_pos_file = op.join(path_to_olga_default_models,
                                    chain_folder,
                                    'V_gene_CDR3_anchors.csv')
        J_anchor_pos_file = op.join(path_to_olga_default_models,
                                    chain_folder,
                                     'J_gene_CDR3_anchors.csv')

            #Load up model based on recomb_type
        #VDJ recomb case --- used for TCRB and IGH
        if recomb_type == 'VDJ':
            genomic_data = load_model.GenomicDataVDJ()
            genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
            generative_model = load_model.GenerativeModelVDJ()
            generative_model.load_and_process_igor_model(marginals_file_name)
            pgen_model = generation_probability.GenerationProbabilityVDJ(generative_model, genomic_data)
            self.genomic_data = genomic_data
            self.generative_model = generative_model
            self.pgen_model = pgen_model

        #VJ recomb case --- used for TCRA and light chain
        elif recomb_type == 'VJ':
            genomic_data = load_model.GenomicDataVJ()
            genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
            generative_model = load_model.GenerativeModelVJ()
            generative_model.load_and_process_igor_model(marginals_file_name)
            pgen_model = generation_probability.GenerationProbabilityVJ(generative_model, genomic_data)
            self.genomic_data = genomic_data
            self.generative_model = generative_model
            self.pgen_model = pgen_model

    def __repr__(self):
        return "tcrdist.pgen.OlgaModel set to recomb_type: '{}' and chain_folder: '{}'".format(self.recomb_type, self.chain_folder)

    
    def compute_aa_cdr3_pgen(self, CDR3_seq, V_usage_mask_in = None, J_usage_mask_in = None):
        """Compute Pgen for the amino acid sequence CDR3_seq.

        ORIGINAL DOC STRING FROM OLGA:

        Conditioned on the V genes/alleles indicated in V_usage_mask_in and the
        J genes/alleles in J_usage_mask_in. (Examples are TCRB sequences/model)

        Parameters
        ----------
        CDR3_seq : str
            CDR3 sequence composed of 'amino acids' -- the standard amino acids,
            plus any custom symbols for an expanded codon alphabet (note the
            standard ambiguous amino acids -- B, J, X, and Z -- are included by
            default).
        V_usage_mask_in : str or list
            An object to indicate which V alleles should be considered. The default
            input is None which returns the list of all productive V alleles.
        J_usage_mask_in : str or list
            An object to indicate which J alleles should be considered. The default
            input is None which returns the list of all productive J alleles.
        print_warnings : bool
            Determines whether warnings are printed or not. Default ON.

        Returns
        -------
        pgen : float
            The generation probability (Pgen) of the sequence

        Examples
        --------
        >>> generation_probability.compute_aa_CDR3_pgen('CAWSVAPDRGGYTF')
        1.5756106696284584e-10
        >>> generation_probability.compute_aa_CDR3_pgen('CAWSVAPDRGGYTF', 'TRBV30*01', 'TRBJ1-2*01')
        1.203646865765782e-10
        >>> generation_probability.compute_aa_CDR3_pgen('CAWXXXXXXXGYTF')
        7.8102586432014974e-05
        """
        pgen_estimate = self.pgen_model.compute_aa_CDR3_pgen(CDR3_seq, V_usage_mask_in, J_usage_mask_in)
        return(pgen_estimate)

    def compute_aa_cdr3_pgens(self, CDR3_seq, V_usage_mask_in = None, J_usage_mask_in = None):
        """
        function for computing many generation probabilities. The assumption is that
        data will come from a pd.DataFrame as lists which can contain either
        strings or None

        Parameters
        ----------

        CDR3_seq : list of strings or None

        Returns
        -------
        pgens : list
            list of floats
        """
        if V_usage_mask_in is None:
            V_usage_mask_in = [None]*len(CDR3_seq)
            if len(CDR3_seq) != len(V_usage_mask_in):
                raise TypeError("len(CDR3_seq) must equal len(V_usage_mask_in)")
        else:
            if len(CDR3_seq) != len(V_usage_mask_in):
                raise TypeError("len(CDR3_seq) must equal len(V_usage_mask_in)")

        if J_usage_mask_in is None:
            J_usage_mask_in = [None]*len(CDR3_seq)
            if len(CDR3_seq) != len(J_usage_mask_in):
                raise TypeError("len(CDR3_seq) must equal len(V_usage_mask_in)")
        else:
            if len(CDR3_seq) != len(J_usage_mask_in):
                raise TypeError("len(CDR3_seq) must equal len(V_usage_mask_in)")

        # Deal with NaN
        l = len(CDR3_seq)
        input_tuples  = [(CDR3_seq[i], V_usage_mask_in[i], J_usage_mask_in[i] ) \
                         for i in range(l)]

        pgen_estimates = [self.compute_aa_cdr3_pgen(x,y,z) for x,y,z in input_tuples]
        return(pgen_estimates)



    def _validate_recomb_type_arg(self):
        ["VDJ","VJ"].index(self.recomb_type)


    def _validate_chain_folder_arg(self):
        """
        checks taht chain argument is valid based on folders in the
        path_to_olga_default_models

        Parameters
        ----------
        chain_folder : string

        Raises
        ------
        ValueError

        """
        valid_chain_folders = [x for x in os.listdir(path_to_olga_default_models)]
        valid_chain_folders_string =  " ".join(map( str, valid_chain_folders))
        try:
            valid_chain_folders.index(self.chain_folder)
        except ValueError:
            raise ValueError('chain_folder arg in OlgaModel must be: {}'.format(valid_chain_folders_string))
        return(1)

    def _validate_chain_folder_with_recomb_type(self):
        """
        check for correct recomb_type cell type match

        TODO: Add DELTA GAMMA and MOUSE RULES
        """
        if self.recomb_type == "VDJ" and self.chain_folder == "human_T_alpha":
            raise ValueError("human_T_alpha default is recomb_type = VJ")
        if self.recomb_type == "VJ" and self.chain_folder == "human_T_beta":
            raise ValueError("human_T_beta default is recomb_type = VDJ")
        return(1)








#
#
#
#
# def olga_pgen(chain_folder = 'human_T_beta',
#               cdr3 = 'CAWSVAPDRGGYTF',
#               v_b = 'TRBV30*01',
#               v_j ='TRBJ1-2*01'):
#     chain_folder = 'human_T_beta'
#     params_file_name = op.join(path_to_olga_default_models,
#                                chain_folder,
#                                'model_params.txt')
#     marginals_file_name = op.join(path_to_olga_default_models,
#                                   chain_folder,
#                                   'model_marginals.txt')
#     V_anchor_pos_file = op.join(path_to_olga_default_models,
#                                 chain_folder,
#                                 'V_gene_CDR3_anchors.csv')
#     J_anchor_pos_file = op.join(path_to_olga_default_models,
#                                 chain_folder,
#                                  'J_gene_CDR3_anchors.csv')
#
#     #Load data
#     genomic_data = load_model.GenomicDataVDJ()
#     genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
#     #Load model
#     generative_model = load_model.GenerativeModelVDJ()
#     generative_model.load_and_process_igor_model(marginals_file_name)
#
#     #Process model/data for pgen computation by instantiating GenerationProbabilityVDJ
#     pgen_model = pgen.GenerationProbabilityVDJ(generative_model, genomic_data)
#
#     #Compute some sequence pgens
#     x = pgen_model.compute_aa_CDR3_pgen('CAWSVAPDRGGYTF', 'TRBV30*01', 'TRBJ1-2*01')
#     assert(np.isclose(x,1.203646865765782e-10))
#     return(x)
