REPORT=sdsspsf
LATEX=pdflatex
BIBTEX=bibtex
RMTEX=rm -f *.dvi *.aux *.log *.blg 

SRCS=$(wildcard *.tex)
REFS=$(wildcard *.bib)
FIGS=$(wildcard FIGURES/*)

all: $(REPORT).pdf

$(REPORT).pdf: $(SRCS) $(REFS) $(FIGS) 
	$(LATEX) $(REPORT)
	$(BIBTEX) $(REPORT)
	$(LATEX) $(REPORT)
	$(LATEX) $(REPORT)
	$(RMTEX)

pdf: $(REPORT).pdf

