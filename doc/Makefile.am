docdir = $(datadir)/doc

docfiles = phyml-manual.tex ./fig/phyrexlog.pdf ./fig/phytimetrace.pdf

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

makeindex: $(MANTEXSRC)
	makeindex $<

phyml-manual.pdf: intro pdflatex makeindex
# phyml-manual.pdf: intro pdflatex
	@echo ""
	@echo "Done."

intro: 
	@echo ""
	@echo ""
	@echo ".: Building documentation :."
	@echo ""
	@echo ""

dist_doc_DATA = $(docfiles)



