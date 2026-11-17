# Copyright 2026 David Lajevec — HEXAPIC
#
# The HEXAPIC website: the home page, the examples and About,
# with the code's reference pages (the hexapic repository's docs/) included
# as docs/ (a copy: `make sync`).  One Sphinx build; the style is the reference pages' own
# (docs/_static/hexapic.css) with the site's additions (_static/site.css).
# Build: README.md.

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if not os.path.isfile(os.path.join(HERE, "docs", "conf.py")):
    raise SystemExit(
        "docs/ is missing: copy the code's docs/ here first, "
        "`make sync HEXAPIC=../hexapic` (README.md)"
    )
sys.path.insert(0, os.path.join(HERE, "docs", "_ext"))

project = "HEXAPIC"
author = "LECAD-PEG, University of Ljubljana"
copyright = "2026, LECAD-PEG, University of Ljubljana"

extensions = [
    "myst_nb",                  # MyST Markdown, and notebooks (not executed at build)
    "sphinx_design",            # cards, grids, tabs, badges, dropdowns
    "sphinx_copybutton",        # a copy button on every code block
    "sphinxcontrib.bibtex",     # citations from docs/references.bib
    "hexapic_site",             # {site}`...` in docs/: here, a link inside this build
    "hexapic_examples",         # examples/gen/: a page per example of the code's examples/ (copied to _examples/)
]
hexapic_examples_src = "_examples"
hexapic_tests_src = "_tests"      # the physics tests' READMEs: rows of the verification's tables
hexapic_site_url = ""           # "": {site} links are internal ones
nb_execution_mode = "off"

bibtex_bibfiles = ["docs/references.bib"]
bibtex_default_style = "unsrt"
bibtex_reference_style = "author_year"

myst_enable_extensions = [
    "amsmath",
    "attrs_inline",
    "colon_fence",
    "deflist",
    "dollarmath",
]
myst_heading_anchors = 3

numfig = True
math_numfig = True
math_eqref_format = "({number})"

exclude_patterns = [
    "_build", ".venv", "README.md", "**/.ipynb_checkpoints",
    "docs/_build", "docs/.venv", "docs/README.md", "docs/_templates",
    "tutorials",                                # not on the site until there are some: its index back in index.md's toctree,
                                                # and here only tutorials/_template.md and tutorials/_drafts
    "_old",                                     # pages replaced by the code's examples, kept until their media is redone
    "_examples", "_tests", "examples/gen/_gallery.md",    # the code's examples (their pages: examples/gen), the gallery (included)
]

html_theme = "pydata_sphinx_theme"
html_title = "HEXAPIC"
html_favicon = "docs/_static/images/favicon.png"
templates_path = ["_templates"]
# the reference pages' static files first (the style, fonts, logos), then the site's
html_static_path = ["docs/_static", "_static"]
html_css_files = ["hexapic.css", "site.css"]
html_js_files = [("hero.js", {"defer": "defer"})]
html_show_sourcelink = False
html_theme_options = {
    "logo": {
        "text": "HEXAPIC",
        "image_light": "docs/_static/images/logo-dark-ink.png",
        "image_dark": "docs/_static/images/logo-light-ink.png",
        "alt_text": "HEXAPIC",
    },
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/LeCAD-PEG/hexapic",
            "icon": "fa-brands fa-github",
        },
    ],
    "navbar_align": "left",
    "header_links_before_dropdown": 6,
    "show_toc_level": 2,
    "footer_start": [],                 # no footer: the site's pages end with their content
    "footer_end": [],
}
