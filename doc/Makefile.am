
docdir = $(datadir)/doc

docfiles = phyml-manual.tex

if HAVE_PDFLATEX
docfiles += phyml-manual.pdf
endif

MANNAME = phyml-manual
MANTEXSRC = $(MANNAME).tex
MANIDX = $(MANNAME).idx

CLEANFILES = $(MANNAME).pdf $(MANNAME).log $(MANNAME).idx $(MANNAME).out \
 $(MANNAME).toc $(MANNAME).aux $(MANNAME).ilg $(MANNAME).ind  $(MANNAME).dvi 

latex: $(MANTEXSRC)
	latex $<

pdflatex: $(MANTEXSRC)
	pdflatex $<

makeindex: $(MANIDX)
	makeindex $<

phyml-manual.pdf: intro latex makeindex pdflatex
	@echo ""
	@echo "Done."

intro: 
	@echo ""
	@echo ""
	@echo ".: Building documentation :."
	@echo ""
	@echo ""

dist_doc_DATA = $(docfiles)



