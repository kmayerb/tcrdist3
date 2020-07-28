"""
Contains the RefGeneSet class

Contains the TCRGene class
"""

from .paths import path_to_db
from collections import OrderedDict
from . import translation
import pandas as pd
import os.path as op

class RefGeneSet:
    """
    RefGeneSet is a class for holding a set of reference TCR sequences

    Attributes
    ----------
    db_file : string
        specifies the name of the file in path_db with reference TCRS
        (e.g alphabeta_db.tsv)
    all_genes : OrderedDict
        contains library of all reference TCRs represented as a dictionary
        of TCRGene instances

    Methods
    -------
    generate_all_genes()
        generates the library of all reference TCRs represented as a dictionary
        of TCRGene instances

    Notes
    -----

    all_genes = RefGeneSet(db_file = "alphabeta_db.tsv").all_genes

    would be equivalent to previously used:

    from all_genes import all_genes


    Information is now callable as follows:
    attribute = 'cdrs'
    RefGeneSet.all_genes['mouse']['TRAV1*01'].__dict__[attribute]
    """
    def __init__(self, db_file):
        self.db_file = db_file
        self.all_genes = self.generate_all_genes(self.db_file)

    def generate_all_genes(self, db_file = "alphabeta_db.tsv"):
        """
        create the reference dataset that will be used throughout tcrdist
        see the TCRGene class below for details

        Parameters
        ----------
        db_file : string
            "alphabeta_db.tsv"

        """
        db_file = op.join(path_to_db, self.db_file)
        assert op.exists(db_file)

        all_genes = OrderedDict()

        db_df = pd.read_csv(db_file, delimiter='\t')

        for rowi, row in db_df.iterrows():
            g = TCRGene( row )
            if g.organism not in all_genes:
                all_genes[g.organism] = OrderedDict() # map from id to TCR_Gene objects
            all_genes[g.organism][g.id] = g
        return(all_genes)


class TCRGene:
    """
    TCRGene is a class for holding information about a single TCR sequence representative

    Parameters
    ----------
    l : list


    Attributes
    ----------
    id : string

    organism : string
        human or mouse
    chain : string
        A or B !!!! this is not the same in TCRrep calls
    region : string
        V or J !!!! this is not the same in TCRrep calls
    protseq: string
        raw protein sequence
    nucseq : string

    alseq : string
        aligned protein sequence
    frame : int

    nucseq_offset : int

    cdr_columns_str : strings
        '28-39;57-66;82-88;106-111',
    cdr_columns : list
        [[28, 39], [57, 66], [82, 88], [106, 111]],
    cdrs : list
        list of cdrs present based on string in the input file l['cdrs']
    cdrs_extracted : list
        list of cdrs extracted based on cdr_column index with gaps retained
    cdrs_aligned : list
        list of cdrs with gaps retained
    cdrs_no_gaps: list
        list of cdrs with gaps removed


    Methods
    -------
    extract_aligned_cdrs()
        wraps _extract_aligned_cdrs
    extract_cdrs_without_gaps()
        wraps _extract_aligned_cdrs

    """
    cdrs_sep = ';'
    gap_character = '.'

    def __init__( self, l):
        self.id       = l['id']
        self.organism = l['organism']
        self.chain    = l['chain']
        self.region   = l['region']
        self.nucseq   = l['nucseq']
        self.alseq    = l['aligned_protseq']
        self.cdr_columns_str = l['cdr_columns']
        self.frame    = l['frame']
        self.rep      = l['id'] # THIS IS PROBABLY NOT ACCUTATE SEE LINE 170 of all genes
        self.mm1_rep  = l['id'] # THIS IS PROBABLY NOT ACCURATE SEE LINES 171 of all genes
        self.nucseq_offset = self.frame - 1 ## 0, 1 or 2 (0-indexed for python)
        self.protseq = translation.get_translation( self.nucseq, self.frame )[0]

        if pd.isnull(l['cdrs']):
            self.cdrs = []
            self.cdr_columns = []
        else:
            self.cdrs = l['cdrs'].split(self.cdrs_sep)
            ## these are still 1-indexed !!!!!!!!!!!!!!
            self.cdr_columns = [ list(map( int, x.split('-'))) for x in l['cdr_columns'].split(self.cdrs_sep) ]


            self.cdrs_extracted = [ self.alseq[ x[0]-1 : x[1] ] for x in self.cdr_columns ]
            self.cdrs_aligned   = self.extract_aligned_cdrs()
            self.cdrs_no_gaps   = self.extract_cdrs_without_gaps()

        # Defensive Assert Statements
        assert self.frame in [1, 2, 3]
        if self.cdrs:
            assert self.cdrs == self.cdrs_extracted # confirm that cdrs in the input file match expectation
            assert self.cdrs_aligned == self.cdrs_extracted

    def extract_aligned_cdrs(self):
        return _extract_aligned_cdrs(self.alseq, self.cdr_columns_str)

    def extract_cdrs_without_gaps(self):
        return _extract_cdrs_without_gaps(self.alseq, self.cdr_columns_str)


def _extract_aligned_cdrs(aligned_protseqs, cdr_columns):
    """

    from db use the aligned_protseq and cdr_column to return the aligned cdrs with gaps preserved


    Parameters
    ----------
    aligned_protseqs : string
        example: 'GQGVEQ.P.AKLMSVEGTFARVNCTYSTSG......FNGLSWYQQREGQAPVFLSYVVL....DGLKDS.....GHFSTFLSRSN.GYSYLLLTELQIKDSASYLCAVR..'
    cdr_columns : string
        example: '28-39;57-66;82-88;106-111')

    Returns
    -------
    cdrs : list


    Examples
    --------
    >>> _extract_aligned_cdrs('GQGVEQ.P.DNLMSVEGTFARVNCTYSTSG......FNGLSWYQQREGHAPVFLSYVVL....DGLKDS.....GHFSTFLSRSN.GYSYLLLTELQIKDSASYLCAVR..',
 '28-39;57-66;82-88;106-111')
    ['TSG......FNG', 'VVL....DGL', 'SRSN.GY', 'CAVR..']
    """
    pos = cdr_columns.split(";")
    pos2 = [list(map(int,x.split("-"))) for x in pos]
    cdrs = [aligned_protseqs[x[0]-1:x[1]] for x in pos2]
    return(cdrs)

def _extract_cdrs_without_gaps(aligned_protseqs, cdr_columns):
    """

    from db use the aligned_protseq and cdr_column to return the aligned cdrs with gaps removed


    Parameters
    ----------
    aligned_protseqs : string
        example: 'GQGVEQ.P.AKLMSVEGTFARVNCTYSTSG......FNGLSWYQQREGQAPVFLSYVVL....DGLKDS.....GHFSTFLSRSN.GYSYLLLTELQIKDSASYLCAVR..'
    cdr_columns : string
        example: '28-39;57-66;82-88;106-111')

    Returns
    -------
    cdrs : list


    Examples
    --------
    >>> _extract_cdrs_without_gaps('GQGVEQ.P.AKLMSVEGTFARVNCTYSTSG......FNGLSWYQQREGQAPVFLSYVVL....DGLKDS.....GHFSTFLSRSN.GYSYLLLTELQIKDSASYLCAVR..',
 '28-39;57-66;82-88;106-111')
    ['TSGFNG', 'VVLDGL', 'SRSNGY', 'CAVR']
    """
    pos = cdr_columns.split(";")
    pos2 = [list(map(int,x.split("-"))) for x in pos]
    cdrs = [aligned_protseqs[x[0]-1:x[1]].replace(".","") for x in pos2]
    return(cdrs)


