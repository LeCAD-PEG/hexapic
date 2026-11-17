# Copyright 2026 David Lajevec — HEXAPIC
# The website's build (README.md).  docs/ is a copy of the code's docs/, made
# by `make sync` (a link will not do: MyST resolves links through it, out of
# the source directory).

HEXAPIC ?= ../hexapic
PY      ?= .venv/bin/python
PORT    ?= 8766
SRC      = $(HEXAPIC)/docs

.PHONY: sync html serve deploy clean

sync:
	@test -f $(SRC)/conf.py || { echo "$(SRC) is not the code's docs/ (set HEXAPIC=)"; exit 1; }
	@rsync -a --delete --exclude .venv/ --exclude _build/ --exclude .ipynb_checkpoints/ $(SRC)/ docs/
	@rsync -a --delete --prune-empty-dirs --include '*/' --include '*.md' --include 'case.py' --include 'figures/*.png' --exclude '*' \
	    $(HEXAPIC)/examples/ _examples/
	@rsync -a --delete --prune-empty-dirs --include '*/' --include 'README.md' --exclude '*' \
	    $(HEXAPIC)/tests/physics/ _tests/

html: sync
	$(PY) -m sphinx -W -b html . _build/html

# rebuilt on every save, here and in the code's docs/ (copied in first); a log
# written here is ignored, or every write of it would start a build
serve: sync
	$(PY) -m sphinx_autobuild . _build/html --host 127.0.0.1 --port $(PORT) \
	    --watch $(SRC) --watch $(HEXAPIC)/examples --pre-build "make -s sync HEXAPIC=$(HEXAPIC)" \
	    --ignore docs --ignore .venv --ignore $(SRC)/_build --ignore $(SRC)/.venv \
	    --ignore "*.log" --ignore "*.swp" --ignore "*~"

# the site as it is here and in the code's docs/ and examples/, published: built strictly into _build/deploy and
# committed to the branch gh-pages (a worktree in _deploy) over what was there, its old/ (the MkDocs site) kept;
# then pushed.  GitHub Pages serves gh-pages at https://lecad-peg.github.io/hexapic/ a minute later.
deploy: sync
	@test -z "$$(git status --porcelain)" || { echo "commit gh-pages-develop first: the deploy names its commit"; exit 1; }
	rm -rf _build/deploy
	$(PY) -m sphinx -W -q -b html . _build/deploy
	rm -rf _build/deploy/.doctrees _build/deploy/.buildinfo
	git fetch -q origin gh-pages
	@test -d _deploy || git worktree add -q _deploy gh-pages
	git -C _deploy checkout -q gh-pages && git -C _deploy merge -q --ff-only origin/gh-pages
	cd _deploy && git ls-files | grep -v '^old/' | grep -v '^\.nojekyll$$' | xargs -r git rm -q --
	cp -r _build/deploy/. _deploy/ && touch _deploy/.nojekyll
	git -C _deploy add -A
	git -C _deploy commit -q -m "Deploy gh-pages-develop $$(git rev-parse --short HEAD)" && git -C _deploy push -q origin gh-pages

clean:
	rm -rf _build
