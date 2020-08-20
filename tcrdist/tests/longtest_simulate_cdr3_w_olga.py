import olga.load_model as load_model
import olga.generation_probability as pgen
import olga.sequence_generation as seq_gen
import pandas as pd


def generate_simulated_beta_seqs(params_file_name = 'tcrdist/default_models/human_T_beta/model_params.txt',
                            marginals_file_name = 'tcrdist/default_models/human_T_beta/model_marginals.txt',
                            V_anchor_pos_file ='tcrdist/default_models/human_T_beta/V_gene_CDR3_anchors.csv',
                            J_anchor_pos_file = 'tcrdist/default_models/human_T_beta/J_gene_CDR3_anchors.csv',
                            output_cols = ['cdr3_b_aa', "v_b_gene",'j_b_gene'],
                            n = 100000):
    #Load data
    genomic_data = load_model.GenomicDataVDJ()
    genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
    #Load model
    generative_model = load_model.GenerativeModelVDJ()
    generative_model.load_and_process_igor_model(marginals_file_name)
    seq_gen_model = seq_gen.SequenceGenerationVDJ(generative_model, genomic_data)

    #Generate some random sequences


    vs=[x[0] for x in genomic_data.__dict__['genV']]
    js=[x[0] for x in genomic_data.__dict__['genJ']]
    vs = {i:k for i,k in enumerate(vs)}
    js = {i:k for i,k in enumerate(js)}

    sim_cdr3 = [seq_gen_model.gen_rnd_prod_CDR3()[1:4] for x in range(n)]
    sim_cdr3_long = [(i,vs[v],js[j]) for i,v,j in sim_cdr3 ]

    df = pd.DataFrame(sim_cdr3_long, columns = output_cols)
    return df

def generate_simulated_alpha_seqs(params_file_name = 'tcrdist/default_models/human_T_alpha/model_params.txt',
                            marginals_file_name = 'tcrdist/default_models/human_T_alpha/model_marginals.txt',
                            V_anchor_pos_file ='tcrdist/default_models/human_T_alpha/V_gene_CDR3_anchors.csv',
                            J_anchor_pos_file = 'tcrdist/default_models/human_T_alpha/J_gene_CDR3_anchors.csv',
                            output_cols = ['cdr3_a_aa', "v_a_gene",'j_a_gene'],
                            n = 100000):
    #Load data
    genomic_data = load_model.GenomicDataVJ()
    genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
    #Load model
    generative_model = load_model.GenerativeModelVJ()
    generative_model.load_and_process_igor_model(marginals_file_name)
    seq_gen_model = seq_gen.SequenceGenerationVJ(generative_model, genomic_data)

    #Generate some random sequences
    vs=[x[0] for x in genomic_data.__dict__['genV']]
    js=[x[0] for x in genomic_data.__dict__['genJ']]
    vs = {i:k for i,k in enumerate(vs)}
    js = {i:k for i,k in enumerate(js)}

    sim_cdr3 = [seq_gen_model.gen_rnd_prod_CDR3()[1:4] for x in range(n)]
    sim_cdr3_long = [(i,vs[v],js[j]) for i,v,j in sim_cdr3 ]

    df = pd.DataFrame(sim_cdr3_long, columns = output_cols)
    return df

if __name__ == "__main__":
    """
    Using Olga See: 
    ---------------
    Zachary Sethna, Yuval Elhanati, Curtis G Callan, Aleksandra M Walczak, Thierry Mora
    `Bioinformatics (2019) <https://doi.org/10.1093/bioinformatics/btz035>`_ 
    OLGA: fast computation of generation probabilities of B- and T-cell receptor amino acid sequences and motifs


    Generate 1000K (1M) CDR3s using default Olga Models
    Human (Alpha/Beta) and Mouse (Beta)

    human_T_alpha_sim1000K.csv
    human_T_beta_sim1000K.csv
    mouse_T_beta_sim1000K.csv
    
    contained in: 
    olga_T_alpha_beta_1000K_simulated_cdr3.zip 
    """
    dfb= generate_simulated_beta_seqs(params_file_name = 'tcrdist/default_models/human_T_beta/model_params.txt',
                                marginals_file_name = 'tcrdist/default_models/human_T_beta/model_marginals.txt',
                                V_anchor_pos_file ='tcrdist/default_models/human_T_beta/V_gene_CDR3_anchors.csv',
                                J_anchor_pos_file = 'tcrdist/default_models/human_T_beta/J_gene_CDR3_anchors.csv',
                                output_cols = ['cdr3_b_aa', "v_b_gene",'j_b_gene'], n = 1000000)
    dfb.to_csv('human_T_beta_sim1000K.csv', index = False)

    dfa = generate_simulated_alpha_seqs(params_file_name = 'tcrdist/default_models/human_T_alpha/model_params.txt',
                                marginals_file_name = 'tcrdist/default_models/human_T_alpha/model_marginals.txt',
                                V_anchor_pos_file ='tcrdist/default_models/human_T_alpha/V_gene_CDR3_anchors.csv',
                                J_anchor_pos_file = 'tcrdist/default_models/human_T_alpha/J_gene_CDR3_anchors.csv',
                                output_cols = ['cdr3_a_aa', "v_a_gene",'j_a_gene'],
                                n = 1000000)

    dfa.to_csv('human_T_alpha_sim1000K.csv', index = False)                         

    dfb= generate_simulated_beta_seqs(params_file_name = 'tcrdist/default_models/mouse_T_beta/model_params.txt',
                                marginals_file_name = 'tcrdist/default_models/mouse_T_beta/model_marginals.txt',
                                V_anchor_pos_file ='tcrdist/default_models/mouse_T_beta/V_gene_CDR3_anchors.csv',
                                J_anchor_pos_file = 'tcrdist/default_models/mouse_T_beta/J_gene_CDR3_anchors.csv',
                                output_cols = ['cdr3_b_aa', "v_b_gene",'j_b_gene'], n = 1000000)
    dfb.to_csv('mouse_T_beta_sim1000K.csv', index = False)
