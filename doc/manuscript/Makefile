LATEX = $(wildcard *.tex)
PDF = $(LATEX:.tex=.pdf)
JUNK = $(wildcard *.log *.aux *.pdf *.blg *.bbl *.nav *.aux *.out *.toc)

all: $(PDF)

%.pdf: %.tex biblio.bib
	pdflatex $*.tex
	if ( grep -q citation $*.aux ) ; then \
		bibtex $* ; \
		pdflatex $*.tex ; \
	fi
	pdflatex $*.tex

clean:
	rm -rf $(PDF) $(JUNK)
