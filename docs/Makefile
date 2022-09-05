# Minimal makefile for Sphinx documentation

# You can set these variables from the command line.
SPHINXOPTS    = -j auto
SPHINXBUILD   = sphinx-build
SOURCEDIR     = source
BUILDDIR      =../../besca-docs
PDFBUILDDIR = /tmp
PDF = ../besca-manual.pdf

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# hack to get the riverplot working. The CDN needs to be at the top of the html file
# Note: only works with GNU sed, not the OS X / macOS / FreeBSD sed!
	sed -i '30i <script src="https://cdn.plot.ly/plotly-2.14.0.min.js"></script>' ../../besca-docs/html/auto_examples/plotting/plot_riverplot.html

latexpdf:
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(PDFBUILDDIR)/latex
	#                                          ^^^
	@echo "Running LaTeX files through pdflatex..."
	make -C $(PDFBUILDDIR)/latex all-pdf
	#         ^^^
	cp $(PDFBUILDDIR)/latex/*.pdf $(PDF)
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	@echo "pdflatex finished; see $(PDF)"

buildandcommithtml: html

	cd $(BUILDDIR)/html; git add . ; git commit -m "rebuilt docs"; git push origin gh-pages
