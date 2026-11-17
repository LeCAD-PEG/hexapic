# The HEXAPIC website

This repository builds <https://lecad-peg.github.io/hexapic/>. It holds the
home page, the examples and the About page. The reference
pages of the code are included as `docs/`.

| This repository (`hexapic-website`) | The code repository (`hexapic/docs/`) |
|---|---|
| The home page: `index.md`, `_templates/hero.html`, `_static/hero.js` | Get started, User guide, Visualization |
| The tutorials, not published yet: `tutorials/` | Implementation, Verification, Performance |
| The examples: `examples/`, `_static/videos/`, `_static/images/examples/` | Development, References |
| About: `about/` | The style: `_static/hexapic.css`, fonts, logos |

The examples' pages come from the code's `examples/`: `make sync` copies
their READMEs, figures and scripts to `_examples/`, and the code's extension
`docs/_ext/hexapic_examples.py` writes a page per example
(`examples/gen/<name>.md`) and the gallery of `examples/index.md`, a section
per group, at the build. This repository adds only the media, by the
example's name: `_static/videos/<name>.mp4` and `.jpg`, and
`_static/images/examples/<name>_thumb.png`. The pages written here before that, and
the recordings of HEXAPIC / VIZ they showed, are kept in `_old/`, which is
not built. The videos, posters and card pictures are made by the code's
`examples/video.py` from a run of each example.

The tutorials are not on the site yet. `tutorials/` holds their index, the
template `_template.md` and drafts, and `conf.py` excludes the directory
from the build. To publish the first ones, add `tutorials/index` to the
toctree of `index.md` and exclude only the template and the drafts.

A video or a page written here needs a commit in this repository only. A
page that describes the code changes together with the code, in the code
repository.

## Building

`docs/` is a copy of the `docs/` directory of the code. It is not part of
this repository, so edit those pages in the code repository. Every make
target first copies the directory from `../hexapic`. To use another
location, pass `HEXAPIC=/path/to/hexapic`.

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
make html
make serve
```

`make html` builds into `_build/html` and treats warnings as errors.
`make serve` serves the site at `http://127.0.0.1:8766` and rebuilds it on
every save, both here and in the `docs/` of the code.

The style comes from `docs/_static/hexapic.css` in the code repository. It
is black and white with one accent colour, `--hx-accent`. The file
`_static/site.css` adds the styles of the home page and the tutorials. To
link from `docs/` to a page of this repository, use
`` {site}`Examples <examples/index>` ``.

## Publishing

The site is published by GitHub Pages from the branch `gh-pages` of
`github.com/LeCAD-PEG/hexapic`, at <https://lecad-peg.github.io/hexapic/>.
This repository is that one's branch `gh-pages-develop`: the site's sources.
Work here, commit to `gh-pages-develop`, and publish now and then:

```bash
make deploy
```

`make deploy` builds the site with `-W` (from the code's `docs/` and
`examples/`, `HEXAPIC=../hexapic` by default), puts the build in place of
what `gh-pages` has, keeps its `old/` (the MkDocs site that was published
before), commits and pushes. `gh-pages` holds only built pages: do not edit
it by hand. The site does not rebuild by itself when the code's pages change:
deploy again.

`_old/` (not published, not committed) keeps the pages and the recordings
that were replaced.
