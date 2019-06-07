import os
import sys

sys.path.insert(0, os.path.abspath('../..'))
extensions = [
'sphinx.ext.autodoc',]

project = 'xMCA'
copyright = '2019, Chenfgei He'
author = 'Chengfei He'
version = '0.1'
release = '0.1'
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
pygments_style = 'sphinx'
html_theme = 'alabaster'
html_static_path = ['_static']