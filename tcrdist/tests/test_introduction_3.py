import pytest 



def test_introduction_3():
    import pandas as pd
    from tcrdist.repertoire import TCRrep

    df = pd.read_csv("dash.csv")
    tr = TCRrep(cell_df = df, 
                organism = 'mouse', 
                chains = ['alpha','beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv',
                compute_distances = False)

    from tcrdist.plotting import plot_pairings, _write_svg
   
    svg_PA  = plot_pairings(tr.clone_df.loc[tr.clone_df.epitope == "PA"], 
                cols = ['v_b_gene', 'j_b_gene','j_a_gene', 'v_a_gene'], 
                count_col='count')
      
    svg_NP = plot_pairings(tr.clone_df.loc[tr.clone_df.epitope == "NP"], 
            cols = ['v_b_gene', 'j_b_gene', 'j_a_gene', 'v_a_gene'], 
            count_col='count')
    
    _write_svg(svg_PA, name = "PA_gene_usage_plot.svg", dest = ".")
    
    _write_svg(svg_NP, name = "NP_gene_usage_plot.svg", dest = ".")

    import fishersapi
    fishersapi.fishers_frame(tr.clone_df.loc[tr.clone_df.epitope == "NP"], 
                             col_pairs=[('v_b_gene', 'j_b_gene'),
                                        ('v_a_gene', 'j_a_gene'),
                                        ('v_a_gene', 'v_b_gene'),
                                        ('j_a_gene', 'j_b_gene')])