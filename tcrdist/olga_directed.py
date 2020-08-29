#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module for Monte Carlo generation of sequences from a V(D)J recomb model.

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


This module defines classes to randomly sample a V(D)J recomb model and
assemble CDR3 sequences defined by the sampled recombination events.



The method gen_rnd_prod_CDR3 will generate the actual CDR3 sequences and
outputs both the nucleotide sequence and the amino acid sequence. It will also
output the V and J gene/allele choices which can be discarded or stored as
the user desires. If a 'sequence read' (i.e. a sequence which mimics the
primer positions from an actual sequencing experiment) is needed, the V and J
identities can be used to determine what the sequence is outside of the CDR3
region.

Example
-------
>>> import olga.load_model as load_model
>>> import olga.sequence_generation as seq_gen
>>>
>>> params_file_name = './models/human_T_beta/model_params.txt'
>>> marginals_file_name = './models/human_T_beta/model_marginals.txt'
>>> V_anchor_pos_file ='./models/human_T_beta/V_gene_CDR3_anchors.csv'
>>> J_anchor_pos_file = './models/human_T_beta/J_gene_CDR3_anchors.csv'
>>>
>>> genomic_data = load_model.GenomicDataVDJ()
>>> genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
>>>
>>> generative_model = load_model.GenerativeModelVDJ()
>>> generative_model.load_and_process_igor_model(marginals_file_name)
>>>
>>> seq_gen_model = seq_gen.SequenceGenerationVDJ(generative_model, genomic_data)
>>>
>>> seq_gen_model.gen_rnd_prod_CDR3()
('TGTGCCAGCAGTGAAAAAAGGCAATGGGAAAGCGGGGAGCTGTTTTTT', 'CASSEKRQWESGELFF', 27, 8)
>>> seq_gen_model.gen_rnd_prod_CDR3()
('TGTGCCAGCAGTTTAGTGGGAAGGGCGGGGCCCTATGGCTACACCTTC', 'CASSLVGRAGPYGYTF', 14, 1)
>>> seq_gen_model.gen_rnd_prod_CDR3()
('TGTGCCAGCTGGACAGGGGGCAACTACGAGCAGTACTTC', 'CASWTGGNYEQYF', 55, 13)

@author: zacharysethna
"""
from __future__ import division
import numpy as np
from olga.utils import nt2aa, calc_steady_state_dist
import warnings

class SequenceGenerationVDJ(object):
    """Class of to generate sequences from a VDJ generative model.

    Attributes
    ----------
    CPV : ndarray
        Cumulative probability distribution of V usage
    given_V_CPdelV : ndarray
        Conditional cumulative distribution of the V deletions given V.

    CPDJ : ndarray
        Joint cumulative probability distribution of the D and J usages.
    given_J_CPdelJ : ndarray
        Conditional cumulative distribution of the J deletions given the J.

    given_D_CPdelDldelDr : ndarray
        Joint cumulative probability distribution of the D deletions given the
        D.

    CPinsVD : ndarray
        Cumulative probability distribution of the VD (N1) insertion sequence
        length.
    CPinsDJ : ndarray
        Cumulative probability distribution of the DJ (N2) insertion sequence
        length.
    CRvd : ndarray
        Cumulative Markov transition matrix for the VD insertion junction.
    CRdj : ndarray
        Cumulative Markov transition matrix for the DJ insertion junction.
    C_first_nt_bias_insVD : ndarray
        (4,) array of the cumulative probability distribution of the indentity
        of the first nucleotide insertion for the VD junction.
    C_first_nt_bias_insDJ : ndarray
        (4,) array of the cumulative probability distribution of the indentity
        of the first nucleotide insertion for the DJ junction.

    num_J_genes : int
        Number of J genes/alleles.
    num_delDr_poss : int
        Number of delDr possibilities.

    cutV_genomic_CDR3_segs : list of strings
        List of the V germline nucleotide sequences, trimmed to begin at the
        CDR3 region (includes the conserved C residue) with the maximum number
        of reverse complementary palindromic insertions appended.
    cutD_genomic_CDR3_segs : list of strings
        List of the D germline nucleotide sequences, with the maximum number
        of reverse complementary palindromic insertions appended to both ends.
    cutJ_genomic_CDR3_segs : list of strings
        List of the J germline nucleotide sequences, trimmed to end at the
        CDR3 region (includes the conserved F or W residue) with the maximum
        number of reverse complementary palindromic insertions appended.

    """

    def __init__(self, generative_model, genomic_data):
        """Initialize SequenceGenerationVDJ

        This intialization computes all of the cumulative probability
        distributions that will be needed for efficient Monte Carlo sequence
        generation out of a GenerativeModelVDJ.

        Parameters
        ----------
        generative_model : GenerativeModelVDJ
            VDJ generative model class containing the model parameters.
        genomic_data : GenomicDataVDJ
            VDJ genomic data class containing the V, D, and J germline
            sequences and info.

        """

        self.CPV = (generative_model.PV/np.sum(generative_model.PV)).cumsum()
        self.CPDJ = (generative_model.PDJ/np.sum(generative_model.PDJ)).flatten().cumsum()
        self.CinsVD = (generative_model.PinsVD/np.sum(generative_model.PinsVD)).cumsum()
        self.CinsDJ = (generative_model.PinsDJ/np.sum(generative_model.PinsDJ)).cumsum()


        for V in range(generative_model.PdelV_given_V.shape[1]):
            if np.sum(generative_model.PdelV_given_V[:, V])> 0:
                generative_model.PdelV_given_V[:, V] = generative_model.PdelV_given_V[:, V]/np.sum(generative_model.PdelV_given_V[:, V])

        self.given_V_CPdelV = generative_model.PdelV_given_V.T.cumsum(axis = 1)

        for J in range(generative_model.PdelJ_given_J.shape[1]):
            if np.sum(generative_model.PdelJ_given_J[:, J])> 0:
                generative_model.PdelJ_given_J[:, J] = generative_model.PdelJ_given_J[:, J]/np.sum(generative_model.PdelJ_given_J[:, J])

        self.given_J_CPdelJ = generative_model.PdelJ_given_J.T.cumsum(axis = 1)

        for D in range(generative_model.PdelDldelDr_given_D.shape[2]):
            if np.sum(generative_model.PdelDldelDr_given_D[:,  :, D]) > 0:
                generative_model.PdelDldelDr_given_D[:, :, D] = generative_model.PdelDldelDr_given_D[:, :, D]/np.sum(generative_model.PdelDldelDr_given_D[:, :, D])

        self.given_D_CPdelDldelDr = np.array([ generative_model.PdelDldelDr_given_D[:, :, i].flatten().cumsum() for i in range(generative_model.PdelDldelDr_given_D.shape[2])])


        self.C_Rvd = generative_model.Rvd.T.cumsum(axis = 1)
        self.C_Rdj = generative_model.Rdj.T.cumsum(axis = 1)

        if generative_model.first_nt_bias_insVD is None:
            first_nt_bias_insVD = calc_steady_state_dist(generative_model.Rvd)
        else:
            first_nt_bias_insVD = generative_model.first_nt_bias_insVD

        if generative_model.first_nt_bias_insDJ is None:
            first_nt_bias_insDJ = calc_steady_state_dist(generative_model.Rdj)
        else:
            first_nt_bias_insDJ = generative_model.first_nt_bias_insDJ



        self.C_first_nt_bias_insVD = first_nt_bias_insVD.cumsum()
        self.C_first_nt_bias_insDJ = first_nt_bias_insDJ.cumsum()

        self.num_J_genes = generative_model.PDJ.shape[1]
        self.num_delDr_poss = generative_model.PdelDldelDr_given_D.shape[1]


        self.cutV_genomic_CDR3_segs = genomic_data.cutV_genomic_CDR3_segs
        self.cutD_genomic_CDR3_segs = genomic_data.cutD_genomic_CDR3_segs
        self.cutJ_genomic_CDR3_segs = genomic_data.cutJ_genomic_CDR3_segs

    def gen_rnd_prod_CDR3(self, V = None, J= None, conserved_J_residues = 'FVW',):
        """Generate a productive CDR3 seq from a Monte Carlo draw of the model.

        Parameters
        ----------
        conserved_J_residues : str, optional
            Conserved amino acid residues defining the CDR3 on the J side (normally
            F, V, and/or W)

        Returns
        -------
        ntseq : str
            Productive CDR3 nucleotide sequence
        aaseq : str
            CDR3 amino acid sequence (aaseq = nt2aa(ntseq))
        V_choice : int
            Index of V allele chosen to generate the CDR3 seq
        J_choice : int
            Index of J allele chosen to generate the CDR3 seq

        """
        
        coding_pass = False
        counter = 0
        while ~coding_pass and counter < 20:
            counter = counter + 1
            #print(counter)
            if V is not None:
                recomb_events = self.choose_directed_recomb_events(V = V, J = J)
            #    print(recomb_events)
            else:
                recomb_events = self.choose_random_recomb_events()
            #print(recomb_events)

            V_seq = self.cutV_genomic_CDR3_segs[recomb_events['V']]
            #print(V_seq)
            #This both checks that the position of the conserved C is
            #identified and that the V isn't fully deleted out of the CDR3
            #region
            if len(V_seq) <= max(recomb_events['delV'], 0):
                continue

            D_seq = self.cutD_genomic_CDR3_segs[recomb_events['D']]
            J_seq = self.cutJ_genomic_CDR3_segs[recomb_events['J']]

            #We check that the D and J aren't deleted more than allowed. Note
            #the generative model really should reflect this structure already
            if len(D_seq) < (recomb_events['delDl'] + recomb_events['delDr']) or len(J_seq) < recomb_events['delJ']:
                continue

            V_seq = V_seq[:len(V_seq) - recomb_events['delV']]
            D_seq = D_seq[recomb_events['delDl']:len(D_seq)-recomb_events['delDr']]
            J_seq = J_seq[recomb_events['delJ']:]

            if (len(V_seq)+ len(D_seq) + len(J_seq) + recomb_events['insVD'] + recomb_events['insDJ']) % 3 != 0:
                continue


            insVD_seq = rnd_ins_seq(recomb_events['insVD'], self.C_Rvd, self.C_first_nt_bias_insVD)
            insDJ_seq = rnd_ins_seq(recomb_events['insDJ'], self.C_Rdj, self.C_first_nt_bias_insDJ)[::-1] #have to reverse the DJ seq

            #Translate to amino acid sequence, see if productive
            ntseq = V_seq + insVD_seq + D_seq + insDJ_seq + J_seq
            aaseq = nt2aa(ntseq)

            if '*' not in aaseq and aaseq[0]=='C' and aaseq[-1] in conserved_J_residues:
                return ntseq, aaseq, recomb_events['V'], recomb_events['J'], recomb_events
        
        #warnings.warn(f"After {counter} attemps no productive CDR3 found from V:{V} and J:{J}, delV likely exceeds V_seq: '{V_seq}' see possible issue with cutV_genomic_CDR3_segs () ")
        return None

    def choose_random_recomb_events(self):
        """Sample the genomic model for VDJ recombination events.

        Returns
        -------
        recomb_events : dict
            Dictionary of the VDJ recombination events. These are
            integers determining gene choice, deletions, and number of insertions.

        Example
        --------
        >>> sequence_generation.choose_random_recomb_events()
        {'D': 0, 'J': 13, 'V': 36, 'delDl': 2, 'delDr': 13, 'delJ': 10, 'delV': 5, 'insDJ': 6, 'insVD': 9}

        """

        recomb_events = {}
        recomb_events['V'] = self.CPV.searchsorted(np.random.random())

        #For 2D arrays make sure to take advantage of a mod expansion to find indicies
        DJ_choice = self.CPDJ.searchsorted(np.random.random())
        recomb_events['D'] = DJ_choice//self.num_J_genes
        recomb_events['J'] = DJ_choice % self.num_J_genes


        #Refer to the correct slices for the dependent distributions
        recomb_events['delV'] = self.given_V_CPdelV[recomb_events['V'], :].searchsorted(np.random.random())

        recomb_events['delJ'] = self.given_J_CPdelJ[recomb_events['J'], :].searchsorted(np.random.random())

        delDldelDr_choice = self.given_D_CPdelDldelDr[recomb_events['D'], :].searchsorted(np.random.random())

        recomb_events['delDl'] = delDldelDr_choice//self.num_delDr_poss
        recomb_events['delDr'] = delDldelDr_choice % self.num_delDr_poss

        recomb_events['insVD'] = self.CinsVD.searchsorted(np.random.random())
        recomb_events['insDJ'] = self.CinsDJ.searchsorted(np.random.random())

        return recomb_events


    def choose_directed_recomb_events(self, V = None, J = None, D = None):
        """Sample the genomic model for VDJ recombination events.

        Returns
        -------
        recomb_events : dict
            Dictionary of the VDJ recombination events. These are
            integers determining gene choice, deletions, and number of insertions.

        Example
        --------
        >>> sequence_generation.choose_random_recomb_events()
        {'D': 0, 'J': 13, 'V': 36, 'delDl': 2, 'delDr': 13, 'delJ': 10, 'delV': 5, 'insDJ': 6, 'insVD': 9}

        """

        recomb_events = {}
        if V is None:
            recomb_events['V'] = self.CPV.searchsorted(np.random.random())
        else:
            assert isinstance(V, int)
            recomb_events['V'] = V
        #For 2D arrays make sure to take advantage of a mod expansion to find indicies
        
        if J is None:
            DJ_choice = self.CPDJ.searchsorted(np.random.random())
            recomb_events['D'] = DJ_choice//self.num_J_genes
            recomb_events['J'] = DJ_choice % self.num_J_genes
        else:
            DJ_choice = self.CPDJ.searchsorted(np.random.random())
            recomb_events['D'] = DJ_choice//self.num_J_genes
            recomb_events['J'] = J 


        #Refer to the correct slices for the dependent distributions
        recomb_events['delV'] = self.given_V_CPdelV[recomb_events['V'], :].searchsorted(np.random.random())

        recomb_events['delJ'] = self.given_J_CPdelJ[recomb_events['J'], :].searchsorted(np.random.random())

        delDldelDr_choice = self.given_D_CPdelDldelDr[recomb_events['D'], :].searchsorted(np.random.random())

        recomb_events['delDl'] = delDldelDr_choice//self.num_delDr_poss
        recomb_events['delDr'] = delDldelDr_choice % self.num_delDr_poss

        recomb_events['insVD'] = self.CinsVD.searchsorted(np.random.random())
        recomb_events['insDJ'] = self.CinsDJ.searchsorted(np.random.random())

        return recomb_events
#%%
class SequenceGenerationVJ(object):
    """Class of to generate sequences from a VJ generative model.

    Attributes
    ----------
    CPVJ : ndarray
        Joint cumulative probability distribution of the V and J usages.
    CPdelV_given_V : ndarray
        Conditional cumulative distribution of the V deletions given the V.
    CPdelJ_given_J : ndarray
        Conditional cumulative distribution of the J deletions given the J.

    CPinsVJ : ndarray
        Cumulative probability distribution of the VJ (N) insertion sequence
        length.
    CRvj : ndarray
        Cumulative Markov transition matrix for the VJ insertion junction.
    C_first_nt_bias_insVJ : ndarray
        (4,) array of the cumulative probability distribution of the indentity
        of the first nucleotide insertion for the VD junction.

    num_J_genes : int
        Number of J genes/alleles

    cutV_genomic_CDR3_segs : list of strings
        List of the V germline nucleotide sequences, trimmed to begin at the
        CDR3 region (includes the conserved C residue) with the maximum number
        of reverse complementary palindromic insertions appended.
    cutJ_genomic_CDR3_segs : list of strings
        List of the J germline nucleotide sequences, trimmed to end at the
        CDR3 region (includes the conserved F or W residue) with the maximum
        number of reverse complementary palindromic insertions appended.

    """

    def __init__(self, generative_model, genomic_data):
        """Initialize SequenceGenerationVJ

        This intialization computes all of the cumulative probability
        distributions that will be needed for efficient Monte Carlo sequence
        generation out of a GenerativeModelVJ.

        Parameters
        ----------
        generative_model : GenerativeModelVJ
            VJ generative model class containing the model parameters.
        genomic_data : GenomicDataVJ
            VJ genomic data class containing the V and J germline
            sequences and info.

        """

        self.CPVJ = (generative_model.PVJ/np.sum(generative_model.PVJ)).flatten().cumsum()
        self.CPinsVJ = (generative_model.PinsVJ/np.sum(generative_model.PinsVJ)).cumsum()

        for V in range(generative_model.PdelV_given_V.shape[1]):
            if np.sum(generative_model.PdelV_given_V[:, V])> 0:
                generative_model.PdelV_given_V[:, V] = generative_model.PdelV_given_V[:, V]/np.sum(generative_model.PdelV_given_V[:, V])

        self.given_V_CPdelV = generative_model.PdelV_given_V.T.cumsum(axis = 1)

        for J in range(generative_model.PdelJ_given_J.shape[1]):
            if np.sum(generative_model.PdelJ_given_J[:, J])> 0:
                generative_model.PdelJ_given_J[:, J] = generative_model.PdelJ_given_J[:, J]/np.sum(generative_model.PdelJ_given_J[:, J])

        self.given_J_CPdelJ = generative_model.PdelJ_given_J.T.cumsum(axis = 1)


        self.C_Rvj = generative_model.Rvj.T.cumsum(axis = 1)

        if generative_model.first_nt_bias_insVJ == None:
            first_nt_bias_insVJ = calc_steady_state_dist(generative_model.Rvj)
        else:
            first_nt_bias_insVJ = generative_model.first_nt_bias_insVJ

        self.C_first_nt_bias_insVJ = first_nt_bias_insVJ.cumsum()

        self.num_J_genes = generative_model.PVJ.shape[1]

        self.cutV_genomic_CDR3_segs = genomic_data.cutV_genomic_CDR3_segs
        self.cutJ_genomic_CDR3_segs = genomic_data.cutJ_genomic_CDR3_segs

    def gen_rnd_prod_CDR3(self, V = None, J = None, conserved_J_residues = 'FVW'):
        """Generate a productive CDR3 seq from a Monte Carlo draw of the model.

        Parameters
        ----------
        conserved_J_residues : str, optional
            Conserved amino acid residues defining the CDR3 on the J side (normally
            F, V, and/or W)

        Returns
        -------
        ntseq : str
            Productive CDR3 nucleotide sequence
        aaseq : str
            CDR3 amino acid sequence (aaseq = nt2aa(ntseq))
        V_choice : int
            Index of V allele chosen to generate the CDR3 seq
        J_choice : int
            Index of J allele chosen to generate the CDR3 seq

        """
        
        coding_pass = False
        counter = 0
        while ~coding_pass and counter < 30:
            counter = counter + 1
            #print(counter)
            if V is not None:
                recomb_events = self.choose_directed_recomb_events(V = V, J = J)
            #    print(recomb_events)
            else:
                recomb_events = self.choose_random_recomb_events()

            V_seq = self.cutV_genomic_CDR3_segs[recomb_events['V']]

            #This both checks that the position of the conserved C is
            #identified and that the V isn't fully deleted out of the CDR3
            #region
            if len(V_seq) <= max(recomb_events['delV'], 0):
                continue
            J_seq = self.cutJ_genomic_CDR3_segs[recomb_events['J']]

            #We check that J isn't deleted more than allowed. Note the
            #generative model really should reflect this structure already
            if len(J_seq) < recomb_events['delJ']:
                continue

            V_seq = V_seq[:len(V_seq) - recomb_events['delV']]
            J_seq = J_seq[recomb_events['delJ']:]

            if (len(V_seq)+len(J_seq) + recomb_events['insVJ']) % 3 != 0:
                continue


            insVJ_seq = rnd_ins_seq(recomb_events['insVJ'], self.C_Rvj, self.C_first_nt_bias_insVJ)

            #Translate to amino acid sequence, see if productive
            ntseq = V_seq + insVJ_seq + J_seq
            aaseq = nt2aa(ntseq)
            #print(aaseq)
            if '*' not in aaseq and aaseq[0]=='C' and aaseq[-1] in conserved_J_residues:
                return ntseq, aaseq, recomb_events['V'], recomb_events['J'], recomb_events

        return None

    def choose_random_recomb_events(self):
        """Sample the genomic model for VDJ recombination events.

        Returns
        -------
        recomb_events : dict
            Dictionary of the VDJ recombination events. These are
            integers determining gene choice, deletions, and number of insertions.

        Example
        --------
        >>> sequence_generation.choose_random_recomb_events()
        {'J': 13, 'V': 36, 'delJ': 10, 'delV': 5, 'insVJ': 3}

        """
        recomb_events = {}

        #For 2D arrays make sure to take advantage of a mod expansion to find indicies
        VJ_choice = self.CPVJ.searchsorted(np.random.random())
        recomb_events['V'] = VJ_choice//self.num_J_genes
        recomb_events['J'] = VJ_choice % self.num_J_genes

        #Refer to the correct slices for the dependent distributions
        recomb_events['delV'] = self.given_V_CPdelV[recomb_events['V'], :].searchsorted(np.random.random())

        recomb_events['delJ'] = self.given_J_CPdelJ[recomb_events['J'], :].searchsorted(np.random.random())
        recomb_events['insVJ'] = self.CPinsVJ.searchsorted(np.random.random())

        return recomb_events

    def choose_directed_recomb_events(self, V = None, J = None):
        """Sample the genomic model for VDJ recombination events.

        Returns
        -------
        recomb_events : dict
            Dictionary of the VDJ recombination events. These are
            integers determining gene choice, deletions, and number of insertions.

        Example
        --------
        >>> sequence_generation.choose_random_recomb_events()
        {'J': 13, 'V': 36, 'delJ': 10, 'delV': 5, 'insVJ': 3}

        """

        recomb_events = {}
        VJ_choice = self.CPVJ.searchsorted(np.random.random())
        #For 2D arrays make sure to take advantage of a mod expansion to find indicies
        if V is None:
            recomb_events['V'] = VJ_choice//self.num_J_genes
        else:
            recomb_events['V'] = V

        if J is None:    
            recomb_events['J'] = VJ_choice % self.num_J_genes
        else:
            recomb_events['J'] = J

        #Refer to the correct slices for the dependent distributions
        recomb_events['delV'] = self.given_V_CPdelV[recomb_events['V'], :].searchsorted(np.random.random())

        recomb_events['delJ'] = self.given_J_CPdelJ[recomb_events['J'], :].searchsorted(np.random.random())
        recomb_events['insVJ'] = self.CPinsVJ.searchsorted(np.random.random())

        return recomb_events



#%% Function to get sequence identity (of a given length)
def rnd_ins_seq(ins_len, C_R, CP_first_nt):
    """Generate a random insertion nucleotide sequence of length ins_len.

    Draws the sequence identity (for a set length) from the distribution
    defined by the dinucleotide markov model of transition matrix R.

    Parameters
    ----------
    ins_len : int
        Length of nucleotide sequence to be inserted.
    C_R : ndarray
        (4, 4) array of the cumulative transition probabilities defined by the
        Markov transition matrix R
    CP_first_nt : ndarray
        (4,) array of the cumulative probabilities for the first inserted
        nucleotide

    Returns
    -------
    seq : str
        Randomly generated insertion sequence of length ins_len.

    Examples
    --------
    >>> rnd_ins_seq(7, CP_generative_model['C_Rvd'], CP_generative_model['C_first_nt_bias_insVD'])
    'GATGGAC'
    >>> rnd_ins_seq(7, CP_generative_model['C_Rvd'], CP_generative_model['C_first_nt_bias_insVD'])
    'ACCCCCG'
    >>> rnd_ins_seq(3, CP_generative_model['C_Rvd'], CP_generative_model['C_first_nt_bias_insVD'])
    'GCC'

    """
    nt2num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    num2nt = 'ACGT'

    if ins_len == 0:
        return ''

    seq = num2nt[CP_first_nt.searchsorted(np.random.random())]
    ins_len += -1

    while ins_len > 0:
        seq += num2nt[C_R[nt2num[seq[-1]], :].searchsorted(np.random.random())]
        ins_len += -1

    return seq