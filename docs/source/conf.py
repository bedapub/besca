# -*- coding: utf-8 -*-
#
# -- Path setup --------------------------------------------------------------

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------

project = 'BESCA'
copyright = '2024, BEDA'
author = 'BEDA'

# The short X.Y version
version = '3.0'
# The full version, including alpha/beta/rc tags
release = '3.0.0'


# -- General configuration ---------------------------------------------------

needs_sphinx = '7.0'

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
    'sphinx_automodapi.automodapi',
    'sphinx.ext.mathjax',
]

templates_path = ['_templates']

# Napoleon
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

# Autosummary
autosummary_generate = True
autodoc_member_order = 'bysource'

# Automodapi
automodapi_toctreedirnm = '../source/besca'



todo_include_todos = True

source_suffix = ['.rst']
master_doc = 'index'
language = 'en'
exclude_patterns = [
    '_build', 'Thumbs.db', '.DS_Store',
    '**.ipynb_checkpoints',
    '**bescape_tutorial.ipynb',
    '**adata_to_eset.ipynb',
]

pygments_style = 'sphinx'
autodoc_docstring_signature = True
autodoc_mock_imports = []
# Don't execute doctest blocks during doc build
doctest_test_doctest_blocks = ''


def skip(app, what, name, obj, skip, options):
    if name == "__init__":
        return False
    return skip


# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 5,
    'includehidden': True,
    'titles_only': False,
}

html_context = dict(
    display_github=True,
    github_user='bedapub',
    github_repo='besca',
    github_version='besca3',
    conf_py_path='/docs/source/',
)

html_static_path = ['_static']
html_extra_path = ['tutorials_html']


def setup(app):
    app.connect("autodoc-skip-member", skip)
    app.add_css_file('css/custom.css')


# -- Options for HTMLHelp output ---------------------------------------------

htmlhelp_basename = 'BESCAdoc'

# -- Intersphinx -------------------------------------------------------------

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'scanpy': ('https://scanpy.readthedocs.io/en/stable/', None),
    'anndata': ('https://anndata.readthedocs.io/en/stable/', None),
}
