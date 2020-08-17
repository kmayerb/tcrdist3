import sys
import os.path as op
import os

path_to_scripts = op.dirname(op.realpath(__file__))
path_to_db = op.join(path_to_scripts, 'db')
path_to_default_models = op.join(path_to_scripts, 'default_models')
path_to_base = op.dirname(path_to_scripts)
