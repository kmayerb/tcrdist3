FROM continuumio/miniconda3:4.8.2

RUN apt-get update && apt-get install -y procps && apt-get install -y nano && apt-get -y install gcc && apt-get -y install unzip && apt-get -y install curl && apt-get -y install wget

RUN pip install git+https://github.com/kmayerb/tcrdist3.git@0.1.7

RUN pip install python-levenshtein==0.12.0
RUN pip install pytest
RUN pip install ipython

RUN python -c "from tcrdist.setup_tests import *; download_and_extract_zip_file('dash.zip', dest = '.')"
