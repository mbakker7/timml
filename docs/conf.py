# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
# import os
# import sys

project = "TimML"
copyright = "2023, Mark Bakker"
author = "Mark Bakker"

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
# sys.path.insert(0, os.path.abspath("."))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "autoapi.extension",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "IPython.sphinxext.ipython_console_highlighting",  # lowercase didn't work
    "numpydoc",
    "myst_nb",
    "sphinx_design",
    "sphinx.ext.autosectionlabel",
    "sphinxcontrib.bibtex",
]

# templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_short_title = "timml"
html_css_files = ["css/custom.css"]
html_show_sphinx = True
html_show_copyright = True
htmlhelp_basename = "timmldoc"  # Output file base name for HTML help builder.
html_use_smartypants = True
html_show_sourcelink = True

html_theme_options = {
    "github_url": "https://github.com/mbakker7/timml",
    "use_edit_page_button": True,
    "header_links_before_dropdown": 7,
    # "icon_links": [
    #     {
    #         "name": "GitHub",  # Label for this link
    #         "url": "https://github.com/mbakker7/timml",  # required
    #         "icon": "fab fa-github-square",
    #         "type": "fontawesome",  # Default is fontawesome
    #     }
    # ],
}

html_context = {
    "github_user": "mbakker7",
    "github_repo": "timml",
    "github_version": "master",
    "doc_path": "docs",
}

# -- Napoleon settings ----------------------------------------------------------------

napoleon_include_init_with_doc = False
napoleon_use_param = True
napoleon_type_aliases = {"ml": "timml.Model"}

# -- Autodoc, autosummary, and autosectionlabel settings ------------------------------

autodoc_typehints = "description"
autodoc_typehints_format = "short"
# autosummary_generate = True
# autoclass_content = "class"
autosectionlabel_prefix_document = True

# -- AutoAPI settings -----------------------------------------------------------------
autoapi_dirs = ["../timml"]
autoapi_root = "05api"

# -- Numpydoc settings ----------------------------------------------------------------

numpydoc_class_members_toctree = True
numpydoc_show_class_members = False

# -- Set intersphinx Directories ------------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/devdocs", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
}

# -- myst_nb options ------------------------------------------------------------------

nb_execution_allow_errors = True  # Allow errors in notebooks, to see the error online
nb_execution_mode = "auto"
nb_merge_streams = True
myst_enable_extensions = ["dollarmath", "amsmath"]
myst_dmath_double_inline = True

# -- bibtex options ------------------------------------------------------------------

# Add some settings for bibtex
bibtex_bibfiles = ["06about/publications.bib"]
bibtex_reference_style = "author_year"
