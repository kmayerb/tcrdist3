#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Define classes for V(D)J generative models and genomic data.

    Copyright (C) 2018 Zachary Sethna

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
    
The classes GenomicDataVDJ and GenomicDataVJ are used to format the germline
genomic sequences/data.

The classes GenerativeModelVDJ and GenerativeModelVJ are used to format V(D)J
recombination model parameters.

Normally the genomic data and model parameters are being read from IGoR 
inference files, and V and J anchor files that have been prepared. For this
purpose methods are included into these classes to read in IGoR files and set
the relevant attributes based on the parameters read in. 

IF YOU DO NOT WANT TO USE IGOR INFERENCE FILES (OR MATCH THEIR SYNTAX), YOU 
WILL NEED WRITE A METHOD THAT WILL PARSE SOME OTHER FILE SYNTAX.

For GenomicDataVDJ and GenomicDataVJ the method:

genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)

will read in the IGoR parameter file params_file_name for the gene/allele
names and germline sequences. The V_anchor_pos_file and J_anchor_pos_file
indicate where in the germline sequences the conserved residue that is used to
define the CDR3 sequence (in general a cysteine 'C' for the V and a 
phenylalanine 'F' or a tryptophan 'W' for the J) as well as whether the V or J
gene/allele is a functional gene or a pseudo gene. It must be noted that the 
gene/allele names from params_file_name MUST EXACTLY MATCH the names that can 
be found in V_anchor_pos_file and J_anchor_pos_file (the exact strings must
match) for the files to load properly. If you are using an IGoR inference file
that is not one of the default models shipped, be sure that the anchor files
match your IGoR parameter file. 


For GenerativeModelVDJ and GenerativeModelVJ:

generative_model.load_and_process_igor_model(marginals_file_name)

reads and processes an IGoR marginals file to construct the probability
distributions of a generative V(D)J model. The assumed model factorization of
the models are the same as in the paper:

VDJ model:
    
P(event) = P(V)P(delV|V)P(D,J)P(delDl,delDr|D)P(delJ|J)P(m_1...m_L)P(n_K...n_1)
where the insertion junction probabilities are a dinucleotide Markov model:
P(m_1...m_L) = PinsVD(L) p_0(m_1) \prod_{i=2}^L Rvd(m_i|m_{i-1})
P(n_K...n_1) = PinsDJ(K) q_0(n_1) \prod_{i=2}^K Rdj(n_i|n_{i-1})

VJ model:
P(event) = P(V,J)P(delV|V)P(delJ|J)P(m_1...m_L)
where the insertion junction is again a dinucleotide Markov model:
P(m_1...m_L) = PinsVJ(L) p_0(m_1) \prod_{i=2}^L Rvj(m_i|m_{i-1})    

The attributes of GenerativeModel (and subsequent processed attributes in 
other classes) attempt to preserve the names from the model factorization as
closely as possible (and thus do not follow the proper Python etiquette of all
lower case with underscores for variable/attribute names).


Some of the functions that directly read in IGoR parameter/marginal files
are adapted from Quentin Marcou (author of IGoR).

@author: zacharysethna

"""
from __future__ import print_function, division
import numpy as np
from olga.utils import cutR_seq, cutL_seq, calc_S_single_gene, calc_S_joint_genes, calc_Sins

#%% GenomicData class definitions
class GenomicData(object):
    """Class used to load genomic data that both VDJ and VJ models have.
    
    The class is the parent of classes GenomicDataVDJ and 
    GenomicDataVJ.

    Attributes
    ----------
    genV : list of lists
        List of genomic V information.
    max_delV_palindrome : int
        Maximum number of reverse complementary insertions for a V germline
        sequence.           
    cutV_genomic_CDR3_segs : list of strings
        List of the V germline nucleotide sequences, trimmed to begin at the
        CDR3 region (includes the conserved C residue) with the maximum number
        of reverse complementary palindromic insertions appended.
        
    genJ : list  of lists
        List of genomic J information.       
    max_delJ_palindrome : int
        Maximum number of reverse complementary insertions for a J germline
        sequence
    cutJ_genomic_CDR3_segs : list of strings
        List of the J germline nucleotide sequences, trimmed to end at the
        CDR3 region (includes the conserved F or W residue) with the maximum 
        number of reverse complementary palindromic insertions appended.
        
    """
    def __init__(self):
        """Initialize GenomicData."""
        
        self.genV = None
        self.genJ = None
        
        self.max_delV_palindrome = None
        self.max_delJ_palindrome = None
        
        self.cutV_genomic_CDR3_segs = None
        self.cutJ_genomic_CDR3_segs = None
        
    def anchor_and_curate_genV_and_genJ(self, V_anchor_pos_file, J_anchor_pos_file):
        """Trim V and J germline sequences to the CDR3 region.
        
        Unproductive sequences have an empty string '' for the CDR3 region
        sequence.
        
        Edits the attributes genV and genJ
        
        Parameters
        ----------
        V_anchor_pos_file_name : str
            File name for the conserved residue (C) locations and functionality 
            of each V genomic sequence.
        J_anchor_pos_file_name : str
            File name for the conserved residue (F/W) locations and 
            functionality of each J genomic sequence.
        
        """
        
        V_anchor_pos = load_genomic_CDR3_anchor_pos_and_functionality(V_anchor_pos_file)
        J_anchor_pos = load_genomic_CDR3_anchor_pos_and_functionality(J_anchor_pos_file)
        
        for V in self.genV:
            try:
                if V_anchor_pos[V[0]][0] > 0 and V_anchor_pos[V[0]][1] != 'X': #Check for functionality
                    V[1] = V[2][V_anchor_pos[V[0]][0]:]
                else:
                    V[1] = ''
            except KeyError:
                V[1] = ''
    
        for J in self.genJ:
            try:
                if J_anchor_pos[J[0]][0] > 0 and J_anchor_pos[J[0]][1] != 'X': #Check for functionality
                    J[1] = J[2][:J_anchor_pos[J[0]][0]+3]
                else:
                    J[1] = ''
            except KeyError:
                J[1] = ''
        
    def generate_cutV_genomic_CDR3_segs(self):
        """Add palindromic inserted nucleotides to germline V sequences.
        
        The maximum number of palindromic insertions are appended to the
        germline V segments so that delV can index directly for number of
        nucleotides to delete from a segment.
        
        Sets the attribute cutV_genomic_CDR3_segs.
        
        """
    
        max_palindrome = self.max_delV_palindrome

        self.cutV_genomic_CDR3_segs = []
        for CDR3_V_seg in [x[1] for x in self.genV]:
            if len(CDR3_V_seg) < max_palindrome:
                self.cutV_genomic_CDR3_segs += [cutR_seq(CDR3_V_seg, 0, len(CDR3_V_seg))]
            else:
                self.cutV_genomic_CDR3_segs += [cutR_seq(CDR3_V_seg, 0, max_palindrome)]
                
    def generate_cutJ_genomic_CDR3_segs(self):
        """Add palindromic inserted nucleotides to germline J sequences.
        
        The maximum number of palindromic insertions are appended to the
        germline J segments so that delJ can index directly for number of
        nucleotides to delete from a segment.
        
        Sets the attribute cutJ_genomic_CDR3_segs.
        
        """
        
        max_palindrome = self.max_delJ_palindrome
        self.cutJ_genomic_CDR3_segs = []
        for CDR3_J_seg in [x[1] for x in self.genJ]:
            if len(CDR3_J_seg) < max_palindrome:
                self.cutJ_genomic_CDR3_segs += [cutL_seq(CDR3_J_seg, 0, len(CDR3_J_seg))]
            else:
                self.cutJ_genomic_CDR3_segs += [cutL_seq(CDR3_J_seg, 0, max_palindrome)]

        
class GenomicDataVDJ(GenomicData):
    """Class used to load genomic data of a VDJ model.

    Attributes
    ----------
    genV : list of lists
        List of genomic V information.
    max_delV_palindrome : int
        Maximum number of reverse complementary insertions for a V germline
        sequence.           
    cutV_genomic_CDR3_segs : list of strings
        List of the V germline nucleotide sequences, trimmed to begin at the
        CDR3 region (includes the conserved C residue) with the maximum number
        of reverse complementary palindromic insertions appended.
        
    genD : list of lists
        List of genomic D information.
    max_delDl_palindrome : int
        Maximum number of reverse complementary insertions for a D germline
        sequence from the left (5').
    max_delDr_palindrome : int
        Maximum number of reverse complementary insertions for a D germline
        sequence from the right (3').
    cutD_genomic_CDR3_segs : list of strings
        List of the D germline nucleotide sequences, with the maximum number
        of reverse complementary palindromic insertions appended to both ends.
        
    genJ : list  of lists
        List of genomic J information.       
    max_delJ_palindrome : int
        Maximum number of reverse complementary insertions for a J germline
        sequence
    cutJ_genomic_CDR3_segs : list of strings
        List of the J germline nucleotide sequences, trimmed to end at the
        CDR3 region (includes the conserved F or W residue) with the maximum 
        number of reverse complementary palindromic insertions appended.
        
    """
    
    def __init__(self, params_file_name = None, V_anchor_pos_file = None, J_anchor_pos_file = None):
        """Initialize GenomicDataVDJ.
        
        Parameters
        ----------
        params_file_name : str or None
            File name for a IGOR parameter file.
        V_anchor_pos_file_name : str or None
            File name for the conserved residue (C) locations and functionality 
            of each V genomic sequence.
        J_anchor_pos_file_name : str or None
            File name for the conserved residue (F/W) locations and 
            functionality of each J genomic sequence.
        
        """
        GenomicData.__init__(self)
        
        self.genD = None
        self.delDl_palindrome = None
        self.delDr_palindrome = None 
        self.cutD_genomic_CDR3_segs = None
        
        if all([params_file_name is not None, V_anchor_pos_file is not None, J_anchor_pos_file is not None]):
            self.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
        
    
    def load_igor_genomic_data(self, params_file_name, V_anchor_pos_file, J_anchor_pos_file):
        """Set attributes by loading in genomic data from IGoR parameter file.
        
        Sets attributes genV,  max_delV_palindrome, cutV_genomic_CDR3_segs, 
        genD, max_delDl_palindrome, max_delDr_palindrome, 
        cutD_genomic_CDR3_segs, genJ, max_delJ_palindrome, and 
        cutJ_genomic_CDR3_segs.
        
        Parameters
        ----------
        params_file_name : str
            File name for a IGOR parameter file.
        V_anchor_pos_file_name : str
            File name for the conserved residue (C) locations and functionality 
            of each V genomic sequence.
        J_anchor_pos_file_name : str
            File name for the conserved residue (F/W) locations and 
            functionality of each J genomic sequence.
        
        """
        
        self.genV = read_igor_V_gene_parameters(params_file_name)
        self.genD = read_igor_D_gene_parameters(params_file_name)
        self.genJ = read_igor_J_gene_parameters(params_file_name)
        
        self.anchor_and_curate_genV_and_genJ(V_anchor_pos_file, J_anchor_pos_file)
 
        self.read_VDJ_palindrome_parameters(params_file_name) #Need palindrome info before generating cut_genomic_CDR3_segs

        self.generate_cutV_genomic_CDR3_segs()
        self.generate_cutD_genomic_CDR3_segs()
        self.generate_cutJ_genomic_CDR3_segs()
                
    def generate_cutD_genomic_CDR3_segs(self):
        """Add palindromic inserted nucleotides to germline V sequences.
        
        The maximum number of palindromic insertions are appended to the
        germline D segments so that delDl and delDr can index directly for number 
        of nucleotides to delete from a segment.
        
        Sets the attribute cutV_genomic_CDR3_segs.
        
        """
        max_palindrome_L = self.max_delDl_palindrome
        max_palindrome_R = self.max_delDr_palindrome

        self.cutD_genomic_CDR3_segs = []
        for CDR3_D_seg in [x[1] for x in self.genD]:
            if len(CDR3_D_seg) < min(max_palindrome_L, max_palindrome_R):
                self.cutD_genomic_CDR3_segs += [cutR_seq(cutL_seq(CDR3_D_seg, 0, len(CDR3_D_seg)), 0, len(CDR3_D_seg))]
            else:
                self.cutD_genomic_CDR3_segs += [cutR_seq(cutL_seq(CDR3_D_seg, 0, max_palindrome_L), 0, max_palindrome_R)]
    
    def read_VDJ_palindrome_parameters(self, params_file_name):
        """Read V, D, and J palindrome parameters from file.
        
        Sets the attributes max_delV_palindrome, max_delDl_palindrome,
        max_delDr_palindrome, and max_delJ_palindrome.
    
        Parameters
        ----------
        params_file_name : str
            File name for an IGoR parameter file of a VDJ generative model.
        
        """
        
        params_file = open(params_file_name, 'r')
        
        
        in_delV = False
        in_delDl = False
        in_delDr = False
        in_delJ = False
        
        
        for line in params_file:
            if line.startswith('#Deletion;V_gene;'):
                in_delV = True
                in_delDl = False
                in_delDr = False
                in_delJ = False
            elif line.startswith('#Deletion;D_gene;Three_prime;'):
                in_delV = False
                in_delDl = False
                in_delDr = True
                in_delJ = False
            elif line.startswith('#Deletion;D_gene;Five_prime;'):
                in_delV = False
                in_delDl = True
                in_delDr = False
                in_delJ = False
            elif line.startswith('#Deletion;J_gene;'):
                in_delV = False
                in_delDl = False
                in_delDr = False
                in_delJ = True
            elif any([in_delV, in_delDl, in_delDr, in_delJ]) and line.startswith('%'):
                if int(line.split(';')[-1]) == 0:
                    if in_delV:
                        self.max_delV_palindrome = np.abs(int(line.lstrip('%').split(';')[0]))
                    elif in_delDl:
                        self.max_delDl_palindrome = np.abs(int(line.lstrip('%').split(';')[0]))
                    elif in_delDr:
                        self.max_delDr_palindrome = np.abs(int(line.lstrip('%').split(';')[0]))
                    elif in_delJ:
                        self.max_delJ_palindrome = np.abs(int(line.lstrip('%').split(';')[0]))
            else:
                in_delV = False
                in_delDl = False
                in_delDr = False
                in_delJ = False
                
                
class GenomicDataVJ(GenomicData):
    """Class used to load genomic data of a VJ model.

    Attributes
    ----------
    genV : list of lists
        List of genomic V information.
    max_delV_palindrome : int
        Maximum number of reverse complementary insertions for a V germline
        sequence.           
    cutV_genomic_CDR3_segs : list of strings
        List of the V germline nucleotide sequences, trimmed to begin at the
        CDR3 region (includes the conserved C residue) with the maximum number
        of reverse complementary palindromic insertions appended.
        
    genJ : list  of lists
        List of genomic J information.       
    max_delJ_palindrome : int
        Maximum number of reverse complementary insertions for a J germline
        sequence
    cutJ_genomic_CDR3_segs : list of strings
        List of the J germline nucleotide sequences, trimmed to end at the
        CDR3 region (includes the conserved F or W residue) with the maximum 
        number of reverse complementary palindromic insertions appended.
        
    """
    def __init__(self, params_file_name = None, V_anchor_pos_file = None, J_anchor_pos_file = None):
        """Initialize GenomicDataVJ.
        
        Parameters
        ----------
        params_file_name : str or None
            File name for a IGOR parameter file.
        V_anchor_pos_file_name : str or None
            File name for the conserved residue (C) locations and functionality 
            of each V genomic sequence.
        J_anchor_pos_file_name : str or None
            File name for the conserved residue (F/W) locations and 
            functionality of each J genomic sequence.
        
        """
        GenomicData.__init__(self)
        if all([params_file_name is not None, V_anchor_pos_file is not None, J_anchor_pos_file is not None]):
            self.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
        
    def load_igor_genomic_data(self, params_file_name, V_anchor_pos_file, J_anchor_pos_file):
        """Set attributes by loading in genomic data from IGoR parameter file.
        
        Sets attributes genV, genJ, max_delV_palindrome, max_delJ_palindrome,
        cutV_genomic_CDR3_segs, and cutJ_genomic_CDR3_segs.
        
        Parameters
        ----------
        params_file_name : str
            File name for a IGOR parameter file.
        V_anchor_pos_file_name : str
            File name for the conserved residue (C) locations and functionality 
            of each V genomic sequence.
        J_anchor_pos_file_name : str
            File name for the conserved residue (F/W) locations and 
            functionality of each J genomic sequence.
        
        """
        
        self.genV = read_igor_V_gene_parameters(params_file_name)
        self.genJ = read_igor_J_gene_parameters(params_file_name)
        
        self.anchor_and_curate_genV_and_genJ(V_anchor_pos_file, J_anchor_pos_file)
        
        self.read_igor_VJ_palindrome_parameters(params_file_name)
        
        self.generate_cutV_genomic_CDR3_segs()
        self.generate_cutJ_genomic_CDR3_segs()

            
    def read_igor_VJ_palindrome_parameters(self, params_file_name):
        """Read V and J palindrome parameters from file.
        
        Sets the attributes max_delV_palindrome and max_delJ_palindrome.
    
        Parameters
        ----------
        params_file_name : str
            File name for an IGoR parameter file of a VJ generative model.
        
        """
        params_file = open(params_file_name, 'r')
        
        
        in_delV = False
        in_delJ = False
        
        
        for line in params_file:
            if line.startswith('#Deletion;V_gene;'):
                in_delV = True
                in_delJ = False
            elif line.startswith('#Deletion;J_gene;'):
                in_delV = False
                in_delJ = True
            elif any([in_delV, in_delJ]) and line.startswith('%'):
                if int(line.split(';')[-1]) == 0:
                    if in_delV:
                        self.max_delV_palindrome = np.abs(int(line.lstrip('%').split(';')[0]))
                    elif in_delJ:
                        self.max_delJ_palindrome = np.abs(int(line.lstrip('%').split(';')[0]))
            else:
                in_delV = False
                in_delJ = False
                
#%% Functions for germline sequence load 
#(assumes IGoR files, and the anchor .csv files)
def load_genomic_CDR3_anchor_pos_and_functionality(anchor_pos_file_name):
    """Read anchor position and functionality from file.

    Parameters
    ----------
    anchor_pos_file_name : str
        File name for the functionality and position of a conserved residue 
        that defines the CDR3 region for each V or J germline sequence.
        
    Returns
    -------
    anchor_pos_and_functionality : dict
        Residue anchor position and functionality for each gene/allele.
    
    """
    
    anchor_pos_and_functionality = {}
    anchor_pos_file = open(anchor_pos_file_name, 'r')
    
    first_line = True
    for line in anchor_pos_file:
        if first_line:
            first_line = False
            continue
        
        split_line = line.split(',')
        split_line = [x.strip() for x in split_line]
        anchor_pos_and_functionality[split_line[0]] = [int(split_line[1]), split_line[2].strip().strip('()')]

    return anchor_pos_and_functionality
                
def read_igor_V_gene_parameters(params_file_name):
    """Load raw genV from file.
    
    genV is a list of genomic V information. Each element is a list of three 
    elements. The first is the name of the V allele, the second is the genomic 
    sequence trimmed to the CDR3 region for productive sequences, and the last 
    is the full germline sequence. For this 'raw genV' the middle element is an
    empty string to be filled in later.

    Parameters
    ----------
    params_file_name : str
        File name for a IGOR parameter file.

    Returns
    -------
    genV : list
        List of genomic V information.
    
    """
    params_file = open(params_file_name, 'r')
    
    V_gene_info = {}

    in_V_gene_sec = False
    for line in params_file:
        if line.startswith('#GeneChoice;V_gene;'):
            in_V_gene_sec = True
        elif in_V_gene_sec:
            if line[0] == '%':
                split_line = line[1:].split(';')
                V_gene_info[split_line[0]] = [split_line[1] , int(split_line[2])]
            else:
                break
    params_file.close()
    
    genV = [[]]*len(V_gene_info.keys())
    
    for V_gene in V_gene_info.keys():
        genV[V_gene_info[V_gene][1]] = [V_gene, '', V_gene_info[V_gene][0]]

    return genV

def read_igor_D_gene_parameters(params_file_name):
    """Load genD from file.
    
    genD is a list of genomic D information. Each element is a list of the name
    of the D allele and the germline sequence.

    Parameters
    ----------
    params_file_name : str
        File name for a IGOR parameter file.

    Returns
    -------
    genD : list
        List of genomic D information.
    
    """
    params_file = open(params_file_name, 'r')
    
    D_gene_info = {}

    in_D_gene_sec = False
    for line in params_file:
        if line.startswith('#GeneChoice;D_gene;'):
            in_D_gene_sec = True
        elif in_D_gene_sec:
            if line[0] == '%':
                split_line = line[1:].split(';')
                D_gene_info[split_line[0]] = [split_line[1] , int(split_line[2])]
            else:
                break
    params_file.close()
    
    genD = [[]]*len(D_gene_info.keys())
    
    for D_gene in D_gene_info.keys():
        genD[D_gene_info[D_gene][1]] = [D_gene, D_gene_info[D_gene][0]]

    return genD

def read_igor_J_gene_parameters(params_file_name):
    """Load raw genJ from file.
    
    genJ is a list of genomic J information. Each element is a list of three 
    elements. The first is the name of the J allele, the second is the genomic 
    sequence trimmed to the CDR3 region for productive sequences, and the last 
    is the full germline sequence. For this 'raw genJ' the middle element is an
    empty string to be filled in later.

    Parameters
    ----------
    params_file_name : str
        File name for a IGOR parameter file.

    Returns
    -------
    genJ : list
        List of genomic J information.
    
    """
    params_file = open(params_file_name, 'r')
    
    J_gene_info = {}

    in_J_gene_sec = False
    for line in params_file:
        if line.startswith('#GeneChoice;J_gene;'):
            in_J_gene_sec = True
        elif in_J_gene_sec:
            if line[0] == '%':
                split_line = line[1:].split(';')
                J_gene_info[split_line[0]] = [split_line[1] , int(split_line[2])]
            else:
                break
    params_file.close()
    
    genJ = [[]]*len(J_gene_info.keys())
    
    for J_gene in J_gene_info.keys():
        genJ[J_gene_info[J_gene][1]] = [J_gene, '', J_gene_info[J_gene][0]]

    return genJ

#%% GenerativeModel Class definitions

class GenerativeModelVDJ(object):
    """Class of a VDJ generative model.
    
    Attributes
    ----------
    PV : ndarray
        Probability distribution of V usage
    PdelV_given_V : ndarray
        Conditional distribution of the V deletions given the V allele,
        i.e. P(delV | V)
    
    PDJ : ndarray
        Joint probability distribution of the D and J usages.
    PdelJ_given_J : ndarray
        Conditional distribution of the J deletions given the J allele,
        i.e. P(delJ | J)
    PdelDldelDr_given_D : ndarray
        Joint probability distribution of the D deletions given the D allele,
        i.e. P(delDl, delDr |D)
    
    PinsVD : ndarray
        Probability distribution of the VD (N1) insertion sequence length        
    PinsDJ : ndarray
        Probability distribution of the DJ (N2) insertion sequence length        
    Rvd : ndarray
        Markov transition matrix for the VD insertion junction.
    Rdj : ndarray
        Markov transition matrix for the DJ insertion junction.        
    first_nt_bias_insVD : ndarray
        (4,) array of the probability distribution of the indentity of the 
        first nucleotide insertion for the VD junction.
    first_nt_bias_insDJ : ndarray
        (4,) array of the probability distribution of the indentity of the 
        first nucleotide insertion for the DJ junction.
    SV : float
        Entropy of PV (in bits).
    SdelV : float
        Conditional entropy of PdelV_givenV (in bits).
    SDJ : float
        Entropy of PDJ (in bits).
    SdelD : float
        Conditional entropy of PdelDldelDr_given_D (in bits).
    SdelJ : float
        Conditional entropy of PdelJ_given_J (in bits)
    SinsVD : float
        Entropy of of PinsVD (in bits)
    S_n1_markov : float
        Entropy of the choice of nucleotides at the N1 (VD) junction (in bits).    
    SinsDJ : float
        Entropy of PinsDJ (in bits)
    S_n2_markov : float
        Entropy of the choice of nucleotides at the N2 (DJ) junction (in bits). 
    Sscenario : float
        Total entropy of recombination event scenarios (in bits).
    
    """
    def __init__(self, marginals_file_name = None):
        """Initialize GenerativeModelVDJ.
        
        Parameters
        ----------
        marginals_file_name : str or None
            File name for a IGoR model marginals file.
            
        """
        
        self.PV = None
        self.PinsVD = None
        self.PinsDJ = None
        self.PdelV_given_V = None
        self.PdelJ_given_J = None
        self.PDJ = None
        self.PdelDldelDr_given_D = None
        self.Rvd = None
        self.Rdj = None
        
        #Normal IGoR inference does not infer these parameters, but allow for
        #the formal model
        self.first_nt_bias_insVD = None
        self.first_nt_bias_insDJ = None
        
        self.SV = None
        self.SdelV = None
        self.SDJ = None
        self.SdelD = None
        self.SdelJ = None
        self.SinsVD = None
        self.S_n1_markov = None
        self.SinsDJ = None
        self.S_n2_markov = None
        
        self.Sscenario = None
        
        if marginals_file_name is not None: 
            self.load_and_process_igor_model(marginals_file_name)
        
    def load_and_process_igor_model(self, marginals_file_name):
        """Set attributes by reading a generative model from IGoR marginal file.
        
        Sets attributes PV, PdelV_given_V, PDJ, PdelJ_given_J, 
        PdelDldelDr_given_D, PinsVD, PinsDJ, Rvd, and Rdj.
        Computes entropies SV, SdelV, SDJ, SdelD, SdelJ, SinsVD, S_n1_markov, 
        SinsDJ, S_n2_markov, and Sscenario.
        
        Parameters
        ----------
        marginals_file_name : str
            File name for a IGoR model marginals file.
        
        """
        
        
        raw_model = read_igor_marginals_txt(marginals_file_name)
        
        self.PV = raw_model[0]['v_choice']
        self.PinsVD = raw_model[0]['vd_ins']
        self.PinsDJ = raw_model[0]['dj_ins']
        self.PdelV_given_V = raw_model[0]['v_3_del'].T
        self.PdelJ_given_J = raw_model[0]['j_5_del'].T
        
        #While this class assumes P(V, D, J) factorizes into P(V)*P(D, J), the B cell model
        #infers allowing for the full correlation. Most of the correlation information is due to
        #chromosomal correlation of alleles (i.e. what chromosome each allele is found on).
        #While this information can be interesting for inference purposes, it is discarded here
        #as generally these models may be use for CDR3s from individuals the models weren't inferred
        #from (and thus the chromosomal correlations are incorrect). This also equates the T and B cell
        #models. To reintroduce the chromosomal correlations use V and J usage masks after inferring the 
        #allele identities on each chromosome.
        
        if raw_model[1]['d_gene'] == ['j_choice', 'd_gene']:
        #Factorized P(V, D, J) = P(V)*P(D, J) --- correct for T cell models
            self.PDJ = np.multiply(raw_model[0]['d_gene'].T, raw_model[0]['j_choice'])
        elif raw_model[1]['d_gene'] == ['v_choice', 'j_choice', 'd_gene']:
        #Full P(V, D, J) for B cells --- need to compute the marginal P(D, J)
            PVJ = np.multiply(raw_model[0]['j_choice'].T, raw_model[0]['v_choice']).T
            PVDJ = np.zeros([raw_model[0]['d_gene'].shape[0], raw_model[0]['d_gene'].shape[2], raw_model[0]['d_gene'].shape[1]])
            for v_in in range(raw_model[0]['d_gene'].shape[0]):
                for j_in in range(raw_model[0]['d_gene'].shape[1]):
                    PVDJ[v_in, :, j_in] = PVJ[v_in, j_in]*raw_model[0]['d_gene'][v_in, j_in, :]
            self.PDJ = np.sum(PVDJ, 0)
        else:
            print('Unrecognized model structure -- need to construct P(D, J)')
            return 0
        
        self.PdelDldelDr_given_D = np.transpose(np.multiply(np.transpose(raw_model[0]['d_3_del'], (2, 0, 1)), raw_model[0]['d_5_del']), (2, 0 , 1))
        Rvd_raw = raw_model[0]['vd_dinucl'].reshape((4, 4)).T
        self.Rvd = np.multiply(Rvd_raw, 1/np.sum(Rvd_raw, axis = 0))
        Rdj_raw = raw_model[0]['dj_dinucl'].reshape((4, 4)).T
        self.Rdj = np.multiply(Rdj_raw, 1/np.sum(Rdj_raw, axis = 0))
        
        
        #Compute Entropies
        SV, SdelV = calc_S_single_gene(self.PV, self.PdelV_given_V)
        SDJ, SdelD, SdelJ = calc_S_joint_genes(self.PDJ, self.PdelDldelDr_given_D, self.PdelJ_given_J)
        
        #Current IGoR models don't have first nt biases -- use default of steady state
        SinsVD, S_n1_markov = calc_Sins(self.PinsVD, self.Rvd)
        SinsDJ, S_n2_markov = calc_Sins(self.PinsDJ, self.Rdj)
        
        self.SV = SV
        self.SdelV = SdelV
        self.SDJ = SDJ
        self.SdelD = SdelD
        self.SdelJ = SdelJ
        self.SinsVD = SinsVD
        self.S_n1_markov = S_n1_markov
        self.SinsDJ = SinsDJ
        self.S_n2_markov = S_n2_markov
        
        self.Sscenario = SV + SdelV + SDJ + SdelD + SdelJ + SinsVD + S_n1_markov + SinsDJ + S_n2_markov
        
class GenerativeModelVJ(object):
    """Class of a VJ generative model.
    
    Attributes
    ----------
    PVJ : ndarray
        Joint probability distribution of the V and J usages.
    PdelV_given_V : ndarray
        Conditional distribution of the V deletions given the V allele,
        i.e. P(delV | V)
    PdelJ_given_J : ndarray
        Conditional distribution of the J deletions given the J allele,
        i.e. P(delJ | J)
        
    PinsVJ : ndarray
        Probability distribution of the VJ (N) insertion sequence length              
    Rvj : ndarray
        Markov transition matrix for the VJ insertion junction.       
    first_nt_bias_insVJ : ndarray
        (4,) array of the probability distribution of the indentity of the 
        first nucleotide insertion for the VD junction.
    SVJ : float
        Entropy of PVJ (in bits).
    SdelV : float
        Conditional entropy of PdelV_givenV (in bits).
    SdelJ : float
        Conditional entropy of PdelJ_given_J (in bits).
    SinsVJ : float
        Entropy of PinsVJ (in bits).
    S_n_markov : float
        Entropy of the choice of nucleotides at the N (VJ) junction (in bits).    
    Sscenario : float
        Total entropy of recombination event scenarios (in bits).
    
    """
    def __init__(self, marginals_file_name = None):
        """Initialize GenerativeModelVJ.
        
        Parameters
        ----------
        marginals_file_name : str or None
            File name for a IGoR model marginals file.
            
        """
        self.PVJ = None
        self.PinsVJ = None
        self.PdelV_given_V = None
        self.PdelJ_given_J = None
        self.Rvj = None
        
        #Normal IGoR inference does not infer these parameters, but allow for
        #the formal model
        self.first_nt_bias_insVJ = None
        
        self.SVJ = None
        self.SdelV = None
        self.SdelJ = None
        self.SinsVJ = None
        self.S_n_markov = None
        
        self.Sscenario = None

        if marginals_file_name is not None: 
            self.load_and_process_igor_model(marginals_file_name)

        
        
    def load_and_process_igor_model(self, marginals_file_name):
        """Set attributes by reading a generative model from IGoR marginal file.
        
        Sets attributes PVJ, PdelV_given_V, PdelJ_given_J, PinsVJ, and Rvj.
        
        Parameters
        ----------
        marginals_file_name : str
            File name for a IGoR model marginals file.
        
        """
        
        raw_model = read_igor_marginals_txt(marginals_file_name)
        
        self.PinsVJ = raw_model[0]['vj_ins']
        self.PdelV_given_V = raw_model[0]['v_3_del'].T
        self.PdelJ_given_J = raw_model[0]['j_5_del'].T
        self.PVJ = np.multiply( raw_model[0]['j_choice'].T, raw_model[0]['v_choice']).T
        Rvj_raw = raw_model[0]['vj_dinucl'].reshape((4, 4)).T
        self.Rvj = np.multiply(Rvj_raw, 1/np.sum(Rvj_raw, axis = 0))
        
        
        #Compute Entropies
        SVJ, SdelV, SdelJ = calc_S_joint_genes(self.PVJ, self.PdelV_given_V, self.PdelJ_given_J)
        SinsVJ, S_n_markov = calc_Sins(self.PinsVJ, self.Rvj)
        
        self.SVJ = SVJ
        self.SdelV = SdelV
        self.SdelJ = SdelJ
        self.SinsVJ = SinsVJ
        self.S_n_markov = S_n_markov
        
        self.Sscenario = SVJ + SdelV + SdelJ + SinsVJ + S_n_markov

#%% Function for reading in IGoR marginal files        
def read_igor_marginals_txt(marginals_file_name , dim_names=False):
    """Load raw IGoR model marginals.
    
    Parameters
    ----------
    marginals_file_name : str
        File name for a IGOR model marginals file.

    Returns
    -------
    model_dict : dict
        Dictionary with model marginals.
    dimension_names_dict : dict
        Dictionary that defines IGoR model dependecies.
    
    """
    with open(marginals_file_name,'r') as file:
        #Model parameters are stored inside a dictionary of ndarrays
        model_dict = {}
        dimension_names_dict = {}
        element_name=""
        first = True
        first_dim_line = False
        element_marginal_array = []
        indices_array = []

        for line in file:
            strip_line = line.rstrip('\n') #Remove end of line character
            if strip_line[0]=='@':
                first_dim_line = True
                if not(first):
                    #Add the previous to the dictionnary
                    model_dict[element_name] = element_marginal_array
                else:
                    first = False
				
                element_name = strip_line[1:]

            if strip_line[0]=='$':
                #define array dimensions
                coma_index = strip_line.find(',')
                dimensions = []

                #Get rid of $Dim[
                previous_coma_index = 4
                while coma_index != -1:
                    dimensions.append(int(strip_line[previous_coma_index+1:coma_index]))
                    previous_coma_index = coma_index
                    coma_index = strip_line.find(',',coma_index+1)
			
                #Add last dimension and get rid of the closing bracket 
                dimensions.append(int(strip_line[previous_coma_index+1:-1]))

                element_marginal_array = np.ndarray(shape=dimensions)

            if strip_line[0]=='#':
                if first_dim_line:
                    dimensions_names = []
                    if len(dimensions) > 1:
                        comma_index = strip_line.find(',')
                        opening_bracket_index = strip_line.find('[')
                        while opening_bracket_index != -1:
                            dimensions_names.append(strip_line[opening_bracket_index+1:comma_index])
                            opening_bracket_index = strip_line.find('[',comma_index) 
                            comma_index = strip_line.find(',',opening_bracket_index)
                    first_dim_line = False
                    dimensions_names.append(element_name)
                    dimension_names_dict[element_name] = dimensions_names
                    
                
                #update indices
                indices_array = []
                if len(dimensions) > 1:
                    comma_index = strip_line.find(',')
                    closing_brack_index = strip_line.find(']')					
                    while closing_brack_index != -1:
                        indices_array.append(int(strip_line[comma_index+1:closing_brack_index]))
                        opening_bracket_index = strip_line.find('[',closing_brack_index) 
                        comma_index = strip_line.find(',',opening_bracket_index)
                        closing_brack_index = strip_line.find(']',closing_brack_index+1)
				

            if strip_line[0]=='%':
                #read doubles
                coma_index = strip_line.find(',')
                marginals_values = []

                #Get rid of the %
                previous_coma_index = 0
                while coma_index != -1:
                    marginals_values.append(float(strip_line[previous_coma_index+1:coma_index]))
                    previous_coma_index = coma_index
                    coma_index = strip_line.find(',',coma_index+1)
			
                #Add last dimension and get rid of the closing bracket 
                marginals_values.append(float(strip_line[previous_coma_index+1:]))
                if len(marginals_values)!=dimensions[-1]:
                    print("problem")
                element_marginal_array[tuple(indices_array)] = marginals_values
        model_dict[element_name] = element_marginal_array				
        
        
    return [model_dict,dimension_names_dict]