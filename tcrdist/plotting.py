from .gene_pairing_plot import plot_pairings
import os

__all__ = ['plot_pairings']

def _write_svg(svg, name, dest = '.'):
    fn = os.path.join(dest, name)
    with open(fn, 'w') as fh:
        fh.write(f"{svg}")

def _write_svgs(svgs, name, dest = '.'):
    """ 
    writes a svg string from multiple svgs to a file dest/name 
    
    Parameters 
    ----------
    svgs : list of strs
        list of svgs strings
    name : str
        filename to output svg composite string
    dest : str
        destination folder
    """
    fn = os.path.join(dest, name)
    with open(fn, 'w') as fh:
        fh.write(f"<body>\n")
        
        for svg in svgs:
            fh.write(f"<div>\n")
            fh.write(f"\t{svg}\n")
            fh.write(f"</div>\n")
        fh.write(f"\n</body>")   
