FROM continuumio/miniconda3:4.8.2

RUN apt-get update && apt-get install -y procps && apt-get install -y nano && apt-get -y install gcc

RUN pip install git+https://github.com/kmayerb/tcrdist3.git@0.1.3

RUN pip install python-levenshtein==0.12.0

RUN python -c "from tcrdist.setup_tests import *; download_and_extract_zip_file('dash.zip', dest = '.')"
