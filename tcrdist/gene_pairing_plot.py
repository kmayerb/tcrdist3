import numpy as np
import pandas as pd
import itertools
#from scipy import stats
#import statsmodels.api as sm
#from sklearn.metrics import adjusted_mutual_info_score

from .html_colors import get_html_colors

__all__ = ['plot_pairings']

greek_alpha = '&#x3b1;'
greek_beta  = '&#x3b2;'
greek_gamma = '&#x3b3;'
greek_delta = '&#x3b4;'

segtype2greek_label = { 'VA':'V'+greek_alpha, 'JA':'J'+greek_alpha,
                        'VB':'V'+greek_beta , 'JB':'J'+greek_beta,
                        'VG':'V'+greek_gamma, 'JG':'J'+greek_gamma,
                        'VD':'V'+greek_delta , 'JD':'J'+greek_delta,
                        'v_a_gene':'V'+greek_alpha, 'j_a_gene':'J'+greek_alpha,
                        'v_b_gene':'V'+greek_beta , 'j_b_gene':'J'+greek_beta,
                      'v_g_gene':'V'+greek_gamma, 'j_g_gene':'J'+greek_gamma,
                        'v_d_gene':'V'+greek_delta , 'j_d_gene':'J'+greek_delta}

def _create_svg(cmds, width, height, background_color=None, use_xlink=False):
    out = ''
    extra = '' if not use_xlink else 'xmlns:xlink="http://www.w3.org/1999/xlink"'
    out += '<svg width="{}" height="{}" xmlns="http://www.w3.org/2000/svg" version="1.1" {} >\n'\
              .format(int(width),int(height),extra)
    if background_color:
        out += _rectangle( (0,0), (width, height), background_color, 'white', 0 )
    out += '\n'.join(cmds) + '\n'
    out += '</svg>\n'
    return out

def _rectangle( upper_left, lower_right, fill, stroke, stroke_width=1, dashed=False ):
    stroke_dasharray_style= ';stroke-dasharray:5,5' if dashed else ''
    return '<rect x="{}" y="{}" height="{}" width="{}" style="fill:{};stroke:{};stroke-width:{}{}" />\n'\
        .format( upper_left[0], upper_left[1], lower_right[1]-upper_left[1], lower_right[0]-upper_left[0],
                 fill, stroke, stroke_width, stroke_dasharray_style )

def _make_text( text, lower_left, fontsize,
               extra_tag = None, font_family = "monospace", color = "black", font_weight = "normal" ):
    assert font_weight in ['normal','bold']
    cmd = '<text x="{:.3f}" y="{:.3f}" font-size="{}" font-weight="{}" font-family="{}" fill="{}" xml:space="preserve">{}</text>\n'\
        .format( lower_left[0], lower_left[1], fontsize, font_weight, font_family, color, text )
    return cmd

def _linear_gradient_cmd(gradient_id_counter, x1, y1, x2, y2, offsets, colors, spread_method="pad" ):
    defs_id_prefix = 'TCRGRAD'
    gradient_id_counter += 1
    id = '{}lingrad{:d}'.format(defs_id_prefix, gradient_id_counter)

    stoplines = ''
    assert len(offsets) == len(colors)

    for offset,color in zip( offsets,colors):
        stoplines += """<stop offset="{:.1f}%"   stop-color="{}" stop-opacity="1"/>""".format( offset, color )

    cmd = """<defs>
                    <linearGradient id="{}"
                    x1="{:.1f}%" y1="{:.1f}%"
                    x2="{:.1f}%" y2="{:.1f}%"
                    spreadMethod="{}">
                                    {}
                    </linearGradient>
            </defs>""".format( id, x1, y1, x2, y2, spread_method, stoplines )

    return gradient_id_counter, id, cmd

def _enrichment_glyph_cmds(center, arrow_length, arrow_width, enrichment, add_rectangle=False, box_color='gold' ):
    cmds = []

    num_heads = int(np.floor(np.abs(np.log(enrichment, 2.) ) ))
    if num_heads<1:
        return cmds

    head_sep = 3. * arrow_width
    head_width = 3. * arrow_width
    head_slant = 1.

    if enrichment>1: ## up arrow, text at the top
        line_p0 = [ center[0], center[1] + arrow_length/2. ]
        line_p1 = [ center[0], center[1] - arrow_length/2. ]
        head_step = 1 * head_sep
    else:
        line_p0 = [ center[0], center[1] - arrow_length/2. ]
        line_p1 = [ center[0], center[1] + arrow_length/2. ]
        head_step = -1 * head_sep

    if add_rectangle:
        ## put a nice rounded rectangle around the outside
        rect_x0 = center[0] - head_width - arrow_width
        rect_y0 = center[1] - arrow_length/2 - arrow_width
        rect_width = 2*(head_width+arrow_width)
        rect_height = arrow_length + 2*arrow_width
        rect_round = 3.0 * arrow_width

        cmds.append( '<rect x="{:.3f}" y="{:.3f}" width="{:.3f}" height="{:.3f}" rx="{:.3f}" ry="{:.3f}" stroke="{}" fill="{}"/>'\
                     .format( rect_x0, rect_y0, rect_width, rect_height, rect_round, rect_round,
                              box_color, box_color ))

    ## make the line
    cmds.append( '<line x1="{:.3f}" y1="{:.3f}" x2="{:.3f}" y2="{:.3f}" stroke="black" stroke-width="{}"/>'\
                 .format( line_p0[0], line_p0[1], line_p1[0], line_p1[1], arrow_width ) )

    ## now make the heads
    for head in range(num_heads):
        for xsign in [1,-1]:
            x1 = line_p0[0]
            x2 = x1 + xsign * head_width
            y1 = line_p1[1] + head * head_step
            y2 = y1 + head_slant * head_step

            cmds.append( '<line x1="{:.3f}" y1="{:.3f}" x2="{:.3f}" y2="{:.3f}" stroke="black" stroke-width="{}"/>'\
                         .format( x1,y1,x2,y2, arrow_width ) )

            if False: ## add shifted white line
                y_shift = arrow_width * head_step / head_sep
                cmds.append( '<line x1="{:.3f}" y1="{:.3f}" x2="{:.3f}" y2="{:.3f}" stroke="white" stroke-width="{}"/>'\
                             .format( x1,y1+y_shift,x2,y2+y_shift, arrow_width ) )
    return cmds

def _computeAssociations(df, cols, count_col='Count'):
    res = []
    col1, col2 = cols
    for val1, val2 in itertools.product(df[col1].unique(), df[col2].unique()):
        OR, pvalue, tab = _testAssociation(df, (col1, val1), (col2, val2), count_col=count_col)
        tot = np.sum(tab)
        res.append({'OR':OR,
                    'pvalue':pvalue,
                    'Col0':col1,
                    'Col1':col2,
                    'Val0':val1,
                    'Val1':val2,
                    '00':tab[0, 0],
                    '01':tab[0, 1],
                    '10':tab[1, 0],
                    '11':tab[1, 1]})
    resDf = pd.DataFrame(res)
    # resDf.loc[:, 'qvalue'] = sm.stats.multipletests(resDf['pvalue'].values, method='fdr_bh')[1]
    # resDf = resDf.sort_values(by='pvalue', ascending=True)
    return resDf

def _testAssociation(df, node1, node2, count_col='Count'):
    """Test if the occurence of nodeA paired with nodeB is more/less common than expected.

    Parameters
    ----------
    nodeX : tuple (column, value)
        Specify the node by its column name and the value.

    Returns
    -------
    OR : float
        Odds-ratio associated with the 2x2 contingency table
    pvalue : float
        P-value associated with the Fisher's exact test that H0: OR = 1"""
    col1, val1 = node1
    col2, val2 = node2

    tmp = df[[col1, col2, count_col]].dropna()
    #print(node1, node2, count_col)
    #print(tmp)

    tab = np.zeros((2, 2))
    tab[0, 0] = (((tmp[col1]!=val1) & (tmp[col2]!=val2)) * tmp[count_col]).sum()
    tab[0, 1] = (((tmp[col1]!=val1) & (tmp[col2]==val2)) * tmp[count_col]).sum()
    tab[1, 0] = (((tmp[col1]==val1) & (tmp[col2]!=val2)) * tmp[count_col]).sum()
    tab[1, 1] = (((tmp[col1]==val1) & (tmp[col2]==val2)) * tmp[count_col]).sum()

    # OR, pvalue = stats.fisher_exact(tab)
    OR, pvalue = np.nan, np.nan
    return OR, pvalue, tab

def plot_pairings(cell_df, cols, count_col=None, use_color_gradients=True, other_frequency_threshold=0.01):
    """Diagram depicts the gene-segment pairing structure of the dataset. The four
    genes (cols) are arrayed left to right. Below each gene-type label (eg "VA")
    is a color-stack showing all the TCR clones and how they break down into the different genes for that gene-type. Each clone
    is devoted a constant vertical height in pixels indicated in the text at the top (N pixels in "Nx y-pixel scale"). The curved
    segments joining neighboring gene-stacks show how the two gene distributions pair up, with the thickness of the segments
    corresponding to the number of clones having those two segments (scaled by the indicated y-pixel scale).

    Column names with format: VA or JB (for Valpha and Jbeta) will have the A, B, G, D replaced with its greek character.

    Parameters
    ----------
    cell_df: pd.DataFrame
        Contains gene segment data, one row per clone, optionally with a counts or frequency column.
    cols : list
        List of columns for displaying frequency and pairings, in order from left to right.
    count_col : str
        Optionally provide a count or frequency column for weighting the rows/clones

    Returns
    -------
    raw_svg : str
        Raw SVG txt that can be written to a file.

    Notes
    -----
    For calling this function

    .. code:: python

      import tcrdist as td
      td.plotting.plot_pairings()

    is equivalent to:

    .. code:: python

      import tcrdist as td
      td.gene_pairing_plot.plot_pairings()


    """
    df = cell_df.copy()

    """Not implemented: enrichment should take into account the whole unbiased repertoire"""
    enrichment_glyphs = False

    if count_col is None:
        df = df.assign(Count=1)
        count_col = 'Count'

    params = dict(min_ami_for_colorscale=0.114, # from shuffling experiments
                    max_ami_for_colorscale=0.5,
                    min_entropy_for_colorscale=0.0,
                    max_entropy_for_colorscale=5.0,
                    min_jsd_for_colorscale=0.02259, ## from background vs background comparisons
                    max_jsd_for_colorscale=0.0,
                    min_gene_frequency_for_labels=0.05,
                    reverse_gradients=False)

    gradient_id_counter = 0

    """Pixel units"""
    left_margin = 50
    right_margin = 50
    top_margin = 50
    bottom_margin = 50
    yspacer = 50
    diagram_height = 600
    svg_height = diagram_height + top_margin + bottom_margin

    flat_band = 50
    middle_band = 400
    slope_weight = 100

    ff = 'sans-serif'

    ypixel_scale = diagram_height / df[count_col].sum()
    svg_width = left_margin + right_margin + (len(cols) - 1) * (flat_band + middle_band) + flat_band

    epitope_fontsize = 60
    midpoint = svg_width / 2
    svg_cmds = []
    '''
    svg_cmds.append( _make_text( '{}'.format(epitope, len(tcrs) ),
                                                  [midpoint-0.5*0.6*epitope_fontsize*len(epitope),
                                                   top_margin+epitope_fontsize-20], epitope_fontsize,
                                                  font_family=ff ) )
    '''

    """Compute marginal frequencies for each value within each column"""
    tot = df[count_col].sum()
    counts = {}
    for c in cols:
        """Assign "Other" to all values below a threshold"""
        tmp = df.groupby(c)[count_col].agg(lambda v: np.sum(v) / tot).sort_values(ascending=False)
        ind = tmp.index[tmp < other_frequency_threshold]
        df.loc[df[c].isin(ind), c] = 'Other ' + c
        counts[c] = df.groupby(c)[count_col].agg(lambda v: np.sum(v)).sort_values(ascending=False)
        #print(c)
        #print(counts[c].sum())

    """For each consecutive pair of columns compute 2 x 2 table of frequencies, OR and Fisher's Exact p-value"""
    res = []
    for i in range(len(cols) - 1):
        res.append(_computeAssociations(df, cols=(cols[i], cols[i+1]), count_col=count_col))
    resDf = pd.concat(res, axis=0, ignore_index=True)
    resDf = resDf.sort_values(by='pvalue').set_index(['Val0', 'Val1'])
    #print(resDf[['00', '01', '10', '11']].sum(axis=1))

    for ii in range(len(cols) - 1):
        top_margin = top_margin
        r0 = cols[ii]
        r1 = cols[ii+1]

        x0 = left_margin + ii*( flat_band + middle_band )

        text = segtype2greek_label.get(r0, r0)
        fontsize = 40.
        xtext = x0 + 0.5*flat_band - 0.5*0.6*fontsize*2
        ytext = top_margin + yspacer - 6
        ## hacking
        ytext -= 6
        if ii==0:
            xtext += 8
        svg_cmds.append( _make_text( text, [ xtext, ytext ], fontsize, font_family=ff ) )
        if ii == (len(cols) - 2): ## add the final column label
            text = segtype2greek_label.get(r1, r1)
            xtext = x0+1.5*flat_band-0.5*0.6*fontsize*2+middle_band
            xtext -= 8
            svg_cmds.append( _make_text( text, [ xtext, ytext ], fontsize, font_family=ff ) )

        r0colors = dict(zip(counts[r0].index, get_html_colors(len(counts[r0]))))
        r1colors = dict(zip(counts[r1].index, get_html_colors(len(counts[r1]))))

        a1pixels = {}
        yleft = yspacer + top_margin
        """Nested loops over the alleles (a0, a1) in the pair of columns, r0, r1"""
        for a0, a0count in counts[r0].iteritems():
            y0_right = yspacer + top_margin
            a0color = r0colors[a0]
            for a1, a1count in counts[r1].iteritems():
                a1color = r1colors[a1]
                vj = (a0, a1)
                a1color = r1colors[a1]

                """Frequency of this allele pairing"""
                band_height = resDf.loc[(a0, a1), '11'] * ypixel_scale
                if band_height > 0:
                    """SVG paths with "stroke-width=0" throw errors for image magick"""
                    """Make a spline"""
                    yright = y0_right + a1pixels.get(a1, 0)

                    #line/spline points
                    j_flat_band = flat_band
                    points = [ (np.floor(x0), yleft + 0.5*band_height ),
                               (x0+flat_band, yleft + 0.5*band_height ),
                               (np.ceil(x0 + flat_band + middle_band), yright + 0.5*band_height ),
                               (np.ceil(x0 + flat_band + middle_band + j_flat_band), yright + 0.5*band_height ) ]


                    path1_cmds = 'M {} {} L {} {} M {} {} C {} {}, {} {}, {} {}'\
                        .format( points[0][0], points[0][1], ## start of v-line
                                 points[1][0], points[1][1], ## end point of v-line
                                 points[1][0], points[1][1],
                                 points[1][0] + slope_weight, points[1][1], ## control for spline start
                                 points[2][0] - slope_weight, points[2][1], ## control for spline end
                                 points[2][0], points[2][1] )

                    if use_color_gradients and a0color != a1color:
                        """NOTE: SVG throws an error if trying to apply a gradient to
                        a perfectly horizontal or vertical line due to bounding box issues"""
                        if np.isclose(points[0][1], points[2][1], atol=1e-3):
                            epsilon = 0.005
                        else:
                            epsilon = 0
                        path1a_cmds = 'M {} {} L {} {}'\
                            .format( points[0][0], points[0][1],  ## start of v-line
                                     points[1][0], points[1][1] ) ## end point of v-line
                        svg_cmds.append( '<path d="{}" stroke="{}" stroke-width="{}" fill="none"/>'\
                                                 .format(path1a_cmds, a0color, band_height ) )
                        ## define the gradient
                        path1b_cmds = 'M {} {} C {} {}, {} {}, {} {}'\
                            .format( points[1][0], points[1][1],
                                     points[1][0] + slope_weight, points[1][1], ## control for spline start
                                     points[2][0] - slope_weight, points[2][1], ## control for spline end
                                     points[2][0], points[2][1] + epsilon )
                        #v_line_rhs_fraction = float(flat_band) / (flat_band + middle_band )
                        offsets = [0, 25.0, 75.0, 100]
                        #offsets = [0, 45.0, 55.0, 100]
                        #offsets = [0, 90.0, 99.0, 100]
                        if params['reverse_gradients']:
                            colors = [a1color, a1color, a0color, a0color]
                        else:
                            colors = [a0color, a0color, a1color, a1color]
                        gradient_id_counter, gradient_id, gradient_cmd = _linear_gradient_cmd(gradient_id_counter, 0, 0, 100, 0, offsets, colors )
                        svg_cmds.append( gradient_cmd )
                        svg_cmds.append( '<path d="{}" stroke="url(#{})" stroke-width="{}" fill="none"/>'\
                                                 .format(path1b_cmds, gradient_id, band_height ) )
                    else:
                        svg_cmds.append( '<path d="{}" stroke="{}" stroke-width="{}" fill="none"/>'\
                                                 .format(path1_cmds,a0color, band_height ) )

                    if ii == (len(cols) - 2): ## add the right-most flat band
                        path2_cmds = 'M {} {} L {} {}'\
                            .format( points[2][0], points[2][1], ## start of j-line
                                     points[3][0], points[3][1] ) ## end of j-line

                        svg_cmds.append( '<path d="{}" stroke="{}" stroke-width="{}" fill="none"/>'\
                                                 .format(path2_cmds, a1color, band_height) )

                yleft += band_height
                a1pixels[a1] = a1pixels.get(a1, 0) + band_height
                y0_right += a1count * ypixel_scale

    """Label the alleles in the left stack (and right stack if ii==2)"""
    fontsize = 20
    fontheight = 0.75 * fontsize
    fontwidth =  0.66 * fontsize

    min_height_for_labels = fontheight + 1
    for ii in range(len(cols) - 1):
        x0 = left_margin + ii*( flat_band + middle_band )
        r0 = cols[ii]
        r1 = cols[ii+1]
        r0colors = dict(zip(counts[r0].index, get_html_colors(len(counts[r0]))))
        r1colors = dict(zip(counts[r1].index, get_html_colors(len(counts[r1]))))
        for jj, (r, rcolors) in enumerate([(r0, r0colors), (r1, r1colors)]):
            if jj == 0 and ii > 0:
                continue

            x = x0 + jj*(flat_band + middle_band)
            ystart = yspacer + top_margin
            for a, acount in counts[r].iteritems():
                if acount*ypixel_scale < min_height_for_labels:
                    #print(acount * ypixel_scale)
                    break

                midpoint =  ystart + acount * ypixel_scale * 0.5
                text = a
                lower_left = [x + 2, midpoint + fontheight / 2]

                """Label in white?"""
                bgcolor = rcolors[a]
                textcolor = 'black' if bgcolor!= 'black' else 'white'

                text_width = fontwidth * len(text)
                lower_left_ha = {'left'  : lower_left,
                                 'right' : [ x+flat_band-text_width, midpoint + fontheight/2 ],
                                 'center': [ x+0.5*flat_band-0.5*text_width, midpoint+fontheight/2 ]}
                if jj==0 and ii==0:
                    """Left-most column"""
                    ha = 'left'
                elif jj==1 and ii==(len(cols) - 2):
                    """Right-most column"""
                    ha = 'right'
                else:
                    ha = 'center'
                # print(text, ii, jj, ha, lower_left_ha[ha])
                svg_cmds.append( _make_text( text, lower_left_ha[ha], fontsize, color=textcolor,
                                                              font_family=ff))
                if enrichment_glyphs:
                    """Add an enrichment glyph"""
                    enrich = resDf['OR']
                    if enrich >=2. or enrich <= 0.5:
                        """Add a glyph"""
                        arrow_length = 1.35 * min_height_for_labels
                        arrow_width = 3.5
                        eg_sep = 14.0
                        if 'A' in r:
                            center = [ lower_left_ha[ha][0] + text_width + eg_sep, midpoint ]
                        else:
                            #print rep
                            assert 'B' in r
                            center = [ lower_left_ha[ha][0] - eg_sep, midpoint ]

                        svg_cmds += _enrichment_glyph_cmds(center, arrow_length, arrow_width, enrich)
                ystart += acount * ypixel_scale

    bg_color = None # 'white'
    rawsvg = _create_svg(svg_cmds,
                         width=svg_width,
                         height=svg_height,
                         background_color=bg_color)
    return rawsvg

if __name__ == '__main__':
    np.random.seed(110820)
    n = 50
    df = pd.DataFrame({'VA':np.random.choice(['TRAV14', 'TRAV12', 'TRAV3', 'TRAV23', 'TRAV11', 'TRAV6'], n),
                       'JA':np.random.choice(['TRAJ4', 'TRAJ2', 'TRAJ3','TRAJ5', 'TRAJ21', 'TRAJ13'], n),
                       'VB':np.random.choice(['TRBV14', 'TRBV12', 'TRBV3', 'TRBV23', 'TRBV11', 'TRBV6'], n),
                       'JB':np.random.choice(['TRBJ4', 'TRBJ2', 'TRBJ3','TRBJ5', 'TRBJ21', 'TRBJ13'], n)})
    df = df.assign(Count=1)
    df.loc[:10, 'Count'] = 10
    svg = plot_pairings(df, ['JA', 'VA', 'VB'], count_col='Count')

    """
    import subprocess
    with open('/home/agartlan/gitrepo/test.svg', 'w') as fh:
        fh.write(svg)
    subprocess.check_call('convert -density 200 ../test.svg ../test.png', shell=True)
    """
