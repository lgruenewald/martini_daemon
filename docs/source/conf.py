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

extensions = ["sphinx.ext.autodoc", "autoapi.extension", "sphinx_rtd_theme", "myst_parser"]

html_theme = "sphinx_rtd_theme"

templates_path = ["_templates"]
exclude_patterns = []
autoapi_dirs = ["../../martini_daemon/"]
# concat class and __init__ documentation
autoapi_python_class_content = "init"
# every class should have its own page
autoapi_own_page_level = "class"
# what to document?
autoapi_options = [
    "members",
    "undoc-members",
    "special-members",
    "show-module-summary",
    "imported-members",
    "private-members",
]


# skip private, but not protected members or special members
def skip_private(app, what, name, obj, skip, options):
    elems = name.split(".")
    if any(elem[:2] == "__" and elem[-2:] != "__" for elem in elems):
        skip = True
    also_skip = {
        "_add_mm_force",
        "_rebuild",
        "_remove_mm_force",
        "_set_default_pbc",
        "_bond_context",
    }
    if elems[-1] in also_skip:
        skip = True
    return skip


def setup(sphinx):
    sphinx.connect("autoapi-skip-member", skip_private)


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_static_path = ["_static"]
