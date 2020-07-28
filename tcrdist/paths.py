import sys
import os.path as op
import os

path_to_scripts = op.dirname(op.realpath(__file__))
path_to_db = op.join(path_to_scripts, 'db')
assert op.isdir( path_to_db )
