FROM quay.io/kmayerb/miniconda3:0.0.2

RUN pip install git+https://github.com/kmayerb/tcrdist3.git@0.1.9

RUN pip install python-levenshtein==0.12.0
RUN pip install pytest
RUN pip install ipython

RUN python -c "from tcrdist.setup_tests import *; download_and_extract_zip_file('dash.zip', dest = '.')"
