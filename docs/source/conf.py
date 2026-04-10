# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "Martini Daemon"
copyright = "2026, "
author = ""

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.autodoc", "autoapi.extension", "sphinx_rtd_theme"]

html_theme = "sphinx_rtd_theme"

templates_path = ["_templates"]
exclude_patterns = []
autoapi_dirs = ["../../martini_daemon/"]
# concat class and __init__ documentation
autoapi_python_class_content = "init"
# every class should have its own page
autoapi_own_page_level = "class"
# what to document? default but private members removed
autoapi_options = [
    "members",
    "undoc-members",
    "special-members",
    "show-module-summary",
    "imported-members",
]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_static_path = ["_static"]
