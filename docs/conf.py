# Configuration file for the Sphinx documentation builder.
# last modified by Chen

import os
import sys

sys.path.insert(0, os.path.abspath("../"))


# ---------------------------- Project information ------------------------------
project = "InPheRNo-ChIP"
copyright = "2024, Chen Su, Amin Emad"
author = "Chen Su, Amin Emad"
# The short X.Y version.
version = "1.2"
# The full version, including alpha/beta/rc tags.
release = "1.2"


# ---------------------------- General configuration -----------------------------
extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.graphviz",
]


templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The master toctree document.
master_doc = "index"

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"


# ---------------------------- Options for HTML output ----------------------------
# html_theme = 'furo'
html_theme = "pydata_sphinx_theme"

html_static_path = ["_static"]
