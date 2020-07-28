import numpy as np

"""
This module contains some basic code from the original tcrdist repo that 
supports .get_translation function call in repertoire_db

    self.protseq = translation.get_translation( self.nucseq, self.frame )[0]
"""

"""
from amino_acids.py module
"""
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', \
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


longer_names={'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

short_to_long = {}

for rsd in list(longer_names.keys()):short_to_long[longer_names[rsd]] = rsd

HP = {'I': 0.73, 'F': 0.61, 'V': 0.54, 'L': 0.53, 'W': 0.37,
      'M': 0.26, 'A': 0.25, 'G': 0.16, 'C': 0.04, 'Y': 0.02,
      'P': -0.07, 'T': -0.18, 'S': -0.26, 'H': -0.40, 'E': -0.62,
      'N': -0.64, 'Q': -0.69, 'D': -0.72, 'K': -1.10, 'R': -1.76}

HP['X'] = sum(HP.values())/20.

GES = {'F': -3.7, 'M': -3.4, 'I': -3.1, 'L': -2.8, 'V': -2.6,
       'C': -2.0, 'W': -1.9, 'A': -1.6, 'T': -1.2, 'G': -1.0,
       'S': -0.6, 'P': 0.2,  'Y': 0.7,  'H': 3.0,  'Q': 4.1,
       'N': 4.8,  'E': 8.2,  'K': 8.8,  'D': 9.2,  'R': 12.3}

GES['X'] = sum(GES.values())/20.

## KD values (Kyte-Doolittle) taken from http://web.expasy.org/protscale/pscale/Hphob.Doolittle.html

KD = {'A': 1.8, 'C': 2.5, 'E': -3.5, 'D': -3.5, 'G': -0.4, 'F': 2.8, 'I': 4.5, 'H': -3.2, 'K': -3.9, 'M': 1.9, 'L': 3.8, 'N': -3.5, 'Q': -3.5, 'P': -1.6, 'S': -0.8, 'R': -4.5, 'T': -0.7, 'W': -0.9, 'V': 4.2, 'Y': -1.3}
assert len(KD) == 20

aa_charge =  {}
for a in amino_acids: aa_charge[a] = 0.0
aa_charge['K'] = 1.0
aa_charge['R'] = 1.0
aa_charge['D'] = -1.0
aa_charge['E'] = -1.0
aa_charge['X'] = 0.0

groups = dict( zip( amino_acids, amino_acids ) )

groups['k'] = '[KR]'
groups['d'] = '[DE]'
groups['n'] = '[NQ]'
groups['s'] = '[ST]'
groups['f'] = '[FYWH]'
groups['a'] = '[AGSP]'
groups['v'] = '[VILM]'

"""
from genetic_code.py module
"""

## http://python.genedrift.org/2007/04/19/genetic-code-part-i/
gencode = {
    'ATA': 'I',    #Isoleucine
    'ATC': 'I',    #Isoleucine
    'ATT': 'I',    # Isoleucine
    'ATG': 'M',    # Methionine
    'ACA': 'T',    # Threonine
    'ACC': 'T',    # Threonine
    'ACG': 'T',    # Threonine
    'ACT': 'T',    # Threonine
    'AAC': 'N',    # Asparagine
    'AAT': 'N',    # Asparagine
    'AAA': 'K',    # Lysine
    'AAG': 'K',    # Lysine
    'AGC': 'S',    # Serine
    'AGT': 'S',    # Serine
    'AGA': 'R',    # Arginine
    'AGG': 'R',    # Arginine
    'CTA': 'L',    # Leucine
    'CTC': 'L',    # Leucine
    'CTG': 'L',    # Leucine
    'CTT': 'L',    # Leucine
    'CCA': 'P',    # Proline
    'CCC': 'P',    # Proline
    'CCG': 'P',    # Proline
    'CCT': 'P',    # Proline
    'CAC': 'H',    # Histidine
    'CAT': 'H',    # Histidine
    'CAA': 'Q',    # Glutamine
    'CAG': 'Q',    # Glutamine
    'CGA': 'R',    # Arginine
    'CGC': 'R',    # Arginine
    'CGG': 'R',    # Arginine
    'CGT': 'R',    # Arginine
    'GTA': 'V',    # Valine
    'GTC': 'V',    # Valine
    'GTG': 'V',    # Valine
    'GTT': 'V',    # Valine
    'GCA': 'A',    # Alanine
    'GCC': 'A',    # Alanine
    'GCG': 'A',    # Alanine
    'GCT': 'A',    # Alanine
    'GAC': 'D',    # Aspartic Acid
    'GAT': 'D',    # Aspartic Acid
    'GAA': 'E',    # Glutamic Acid
    'GAG': 'E',    # Glutamic Acid
    'GGA': 'G',    # Glycine
    'GGC': 'G',    # Glycine
    'GGG': 'G',    # Glycine
    'GGT': 'G',    # Glycine
    'TCA': 'S',    # Serine
    'TCC': 'S',    # Serine
    'TCG': 'S',    # Serine
    'TCT': 'S',    # Serine
    'TTC': 'F',    # Phenylalanine
    'TTT': 'F',    # Phenylalanine
    'TTA': 'L',    # Leucine
    'TTG': 'L',    # Leucine
    'TAC': 'Y',    # Tyrosine
    'TAT': 'Y',    # Tyrosine
    'TAA': '_',    # Stop
    'TAG': '_',    # Stop
    'TGC': 'C',    # Cysteine
    'TGT': 'C',    # Cysteine
    'TGA': '_',    # Stop
    'TGG': 'W',    # Tryptophan
    }


## http://www.cmbi.kun.nl/pythoncourse/spy/index.spy?site=python&action=Grand%20Finale&flag=chap

standard = { 'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
             'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
             'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': 'W',
             'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
             'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
             'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
             'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
             'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
             'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
             'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
             'ATA': 'M', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
             'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
             'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
             'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
             'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
             'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}

## http://www.pasteur.fr/formation/infobio/python/ch15.html

code = {'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
        'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
        'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
        'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
        'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
        'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
        'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
        'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
        'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
        'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
        'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
        'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
        'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
        'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
        'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
        'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'
        }

genetic_code = {}
reverse_genetic_code = {}

for codon in gencode:
    lowcodon = codon.lower()
    assert code[ lowcodon ] == gencode[ codon ] or ( gencode[ codon ] == '_' and code[ lowcodon ] == '*' )
    aa = code[ lowcodon ]
    genetic_code[ lowcodon ] = aa
    if aa not in reverse_genetic_code: reverse_genetic_code[aa] = []
    reverse_genetic_code[aa].append( lowcodon )


"""
From translation.py module
"""
base_partner = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n',
                'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                'R': 'Y', 'Y': 'R',
                'S': 'S', 'W': 'W',
                'K': 'M', 'M': 'K',
                '.': '.' }

nucleotide_classes_lower_case = { 'a':'a',
                                  'c':'c',
                                  'g':'g',
                                  't':'t',
                                  'w':'at',
                                  's':'cg',
                                  'k':'gt',
                                  'm':'ac',
                                  'y':'ct',
                                  'r':'ag',
                                  'n':'acgt' }

def reverse_complement( seq ):
    newseq = ''
    L = len(seq)
    for pos in range( L-1, -1, -1 ):
        newseq += base_partner[ seq[ pos ] ]
    assert len( newseq ) == L
    return newseq

def modify_genetic_code(genetic_code):
    bases_plus = list(nucleotide_classes_lower_case.keys())

    for a in bases_plus:
        for b in bases_plus:
            for c in bases_plus:
                codon = a+b+c
                if codon in genetic_code: continue

                aas = []
                for a1 in nucleotide_classes_lower_case[a]:
                    for b1 in nucleotide_classes_lower_case[b]:
                        for c1 in nucleotide_classes_lower_case[c]:
                            aas.append( genetic_code[ a1+b1+c1 ] )
                if min(aas) == max(aas):
                    genetic_code[codon] = aas[0]
                else:
                    genetic_code[codon] = 'X'
    return genetic_code

def get_translation(seq, frame):
    assert frame in [-1, -2, -3, 1, 2, 3, '+1', '+2', '+3', '1', '2', '3', '-1', '-2', '-3']
    frame = int(frame)
    if frame < 0:
        seq = reverse_complement( seq )
    offset = np.abs(frame) - 1
    assert offset in range(3)
    seq = seq[offset:].lower()
    naa = len(seq)//3
    protseq = ''
    codons = []
    for i in range(int(naa)):
        codon = seq[3*i:3*i+3]
        codons.append( codon )
        if '#' in codon:
            protseq += '#'
        else:
            protseq += genetic_code.get( codon, 'X' )
    return protseq, codons

genetic_code = modify_genetic_code(genetic_code)

