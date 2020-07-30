# -*- coding: utf-8 -*-
#
# -- Path setup --------------------------------------------------------------

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------

project = 'BESCA'
copyright = '2020, BEDA'
author = 'BEDA'

# The short X.Y version
version = '2.1'
# The full version, including alpha/beta/rc tags
release = '2.1'


# -- General configuration ---------------------------------------------------

#necessary minimal sphinx version
needs_sphinx = '1.8'

#list of extensions
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary', 
    'sphinx.ext.intersphinx',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.autosectionlabel',
    'sphinx_gallery.gen_gallery',
    'sphinx_automodapi.automodapi',
    'matplotlib.sphinxext.plot_directive',
    'nbsphinx',
    'sphinx.ext.mathjax'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

#configure napoleon extension
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

#autosummary configurations
autosummary_generate = True
autodoc_member_order = 'bysource'

#configure automodapi
automodapi_toctreedirnm = '../source/besca'

#configure sphinx gallery
sphinx_gallery_conf = {
    # path to your examples scripts
    'examples_dirs': '../../besca/examples/gallery_examples',
    # path where to save gallery generated examples
    'gallery_dirs': 'auto_examples',
    'download_all_examples': False,
    'thumbnail_size': (400,400)
}

#include todos
todo_include_todos= True

# The suffix(es) of source filenames.
source_suffix = ['.rst', '.ipynb']

# The master toctree document.
master_doc = 'index'


language = None
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints', '**bescape_tutorial.ipynb', '**adata_to_eset.ipynb']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

autodoc_docstring_signature = True

def skip(app, what, name, obj, skip, options):
    if name == "__init__":
        return False
    return skip

def setup(app):
    app.connect("autodoc-skip-member", skip)

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 5,
    'includehidden': True,
    'titles_only': False
}

html_context = dict(
    display_github=False,      # Integrate GitHub
    github_user='BEDA',       # Username
    github_repo='besca',      # Repo name
    github_version='master',  # Version
    conf_py_path='/docs/',    # Path in the checkout to the docs root
)

html_static_path = ['_static']

html_extra_path = ['tutorials_html']

def setup(app):
    app.add_stylesheet('css/custom.css')

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'BESCAdoc'

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'https://docs.python.org/': None}

# -- Strip output ----------------------------------------------
import nbclean, glob

for filename in glob.glob('**/*.ipynb', recursive=True):
    ntbk = nbclean.NotebookCleaner(filename)
    ntbk.clear('stderr')
    ntbk.save(filename)
