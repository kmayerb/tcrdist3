import sys
import os

path_to_scripts = os.path.dirname(os.path.realpath(__file__))
path_to_db = os.path.join(path_to_scripts, 'db')
path_to_default_models = os.path.join(path_to_scripts, 'default_models')
path_to_base = os.path.dirname(path_to_scripts)
path_to_data = os.path.join(path_to_scripts, 'data')
