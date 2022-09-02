# -*- coding: utf-8 -*-
#
# -- Path setup --------------------------------------------------------------

import os
import sys
from datetime import datetime
import besca

sys.path.insert(0, os.path.abspath("../.."))

# -- Project information -----------------------------------------------------

project = "BESCA"
copyright = f"{datetime.now():%Y}, Roche Pharma Research and Early Development, \
    Pharmaceutical Sciences, Roche Innovation Center Basel, Basel, Switzerland"
author = "PMDA"

# The short X.Y version
version = "2.5"
# The full version, including alpha/beta/rc tags
release = besca.__version__.replace(".dirty", "")


# -- General configuration ---------------------------------------------------

# necessary minimal sphinx version
needs_sphinx = "5.1"
# Warn about broken links
nitpicky = True


# list of extensions
extensions = [
    "myst_parser",
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    "sphinx.ext.todo",
    "sphinx.ext.autosectionlabel",
    "sphinx_gallery.gen_gallery",
    # "sphinx_gallery.load_style",
    "sphinx_automodapi.automodapi",
    "sphinx_copybutton",
    "matplotlib.sphinxext.plot_directive",
    "nbsphinx",
    "sphinx.ext.mathjax",
]

language = "en"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# configure napoleon extension
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

# autosummary configurations
autosummary_generate = True
autodoc_member_order = "bysource"

typehints_defaults = "braces"

# configure automodapi
automodapi_toctreedirnm = "../source/besca"

# configure sphinx gallery
sphinx_gallery_conf = {
    # path to your examples scripts
    "examples_dirs": "../../besca/examples/gallery_examples",
    # path where to save gallery generated examples
    "gallery_dirs": "auto_examples",
    "download_all_examples": False,
    "thumbnail_size": (400, 400),
}

# include todos
todo_include_todos = True

# The suffix(es) of source filenames.
# note: do not add .ipynb when nbsphinx is enabled, otherwise you get the "missing title" error
# https://stackoverflow.com/questions/55297443/including-notebook-with-nbsphinx-fails
source_suffix = [".rst", ".txt" ".md"]

nbsphinx_allow_errors = False

# The master toctree document.
master_doc = "index"


exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**.ipynb_checkpoints",
    "**bescape_tutorial.ipynb",
    "**adata_to_eset.ipynb",
]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

autodoc_docstring_signature = True


def skip(app, what, name, obj, skip, options):
    if name == "__init__":
        return False
    return skip


def setup(app):
    app.connect("autodoc-skip-member", skip)
    app.add_css_file("css/custom.css")
    app.warningiserror = False


# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"

html_theme_options = {
    "logo_only": False,
    "display_version": True,
    "prev_next_buttons_location": "bottom",
    "style_external_links": False,
    # Toc options
    "collapse_navigation": True,
    "sticky_navigation": True,
    "navigation_depth": 5,
    "includehidden": True,
    "titles_only": False,
}

html_context = dict(
    display_github=False,  # Integrate GitHub
    github_user="BEDA",  # Username
    github_repo="besca",  # Repo name
    github_version="master",  # Version
    conf_py_path="/docs/",  # Path in the checkout to the docs root
)

html_static_path = ["_static"]

html_extra_path = ["tutorials_html"]


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = "BESCAdoc"

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {"https://docs.python.org/": None}

# -- Strip output ----------------------------------------------
import nbclean
import glob

for filename in glob.glob("**/*.ipynb", recursive=True):
    ntbk = nbclean.NotebookCleaner(filename)
    ntbk.clear("stderr")
    ntbk.save(filename)
