# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

from functools import partial

from docutils import nodes
from sphinx.application import Sphinx

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "MORESCA"
copyright = "2025, Matthias Bruhns, Jan T. Schleicher"
author = "Matthias Bruhns, Jan T. Schleicher"
release = "0.3.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["myst_parser", "sphinx_copybutton"]

templates_path = ["_templates"]
exclude_patterns = []

source_suffix = {".rst": "restructuredtext", ".txt": "markdown", ".md": "markdown"}

suppress_warnings = ["myst.header"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
html_theme_options = dict(
    use_repository_button=True,
    repository_url="https://github.com/claassenlab/MORESCA/",
    repository_branch="main",
)
html_title = "MORESCA"


def setup(app: Sphinx) -> None:
    """App setup hook."""
    app.add_generic_role("small", partial(nodes.inline, classes=["small"]))
