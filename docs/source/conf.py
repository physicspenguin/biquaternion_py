# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "biquaternion_py"
copyright = "2023, Daren Thimm"
author = "Daren Thimm"
release = "0.1.2"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx.ext.duration",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "numpydoc",
    "nbsphinx",
]

templates_path = ["_templates"]
exclude_patterns = []

# Numpydoc config
# exclude inherited class members
numpydoc_show_inherited_class_members = False
numpydoc_show_class_members = False

# viewcode config
viewcode_line_numbers = True

# nbsphinx config
hightlight_language = "ipython"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
