###############################################################################
# draco/doc/autodoc Makefile
# John M. McGhee
# Mon Aug 9  18:07:47 1999
# Time-stamp: <99/08/09 18:10:19 mcghee>
###############################################################################
# @> Makefile for Draco autodoc document
###############################################################################

DRACO := ../..
DRACOTEXINPUTS := ${DRACO}/doc/tex
TEXINPUTSENVS := TEXINPUTS="${DRACOTEXINPUTS}:${TEXINPUTS}"

LATEXENVS := ${TEXINPUTSENVS}
LATEX := ${LATEXENVS} latex

BIBTEXENVS := BSTINPUTS="${DRACOTEXINPUTS}:${BSTINPUTS}"
BIBTEX := ${BIBTEXENVS} bibtex

MAKEIDX := makeindex

.PHONY : all install clean dvi postscript clean-tex index bib tex

all: postscript

install: all
	@echo install not ready

postscript: draco-autodoc.ps

dvi: draco-autodoc.dvi

draco-autodoc.ps: draco-autodoc.dvi
	dvips -o draco-autodoc.ps draco-autodoc.dvi

draco-autodoc.dvi: draco-autodoc.tex 
	${LATEX} draco-autodoc
	${BIBTEX} draco-autodoc
	${LATEX} draco-autodoc
	${LATEX} draco-autodoc

index: 
	${MAKEIDX} draco-autodoc

bib:
	${BIBTEX} draco-autodoc

tex:
	${LATEX} draco-autodoc

clean:
	rm -f draco-autodoc.dvi draco-autodoc.ps draco-autodoc.log draco-autodoc.aux
	rm -f draco-autodoc.blg draco-autodoc.bbl
	rm -f draco-autodoc.idx draco-autodoc.ilg draco-autodoc.ind
	rm -f draco-autodoc.toc draco-autodoc.lof draco-autodoc.lot

clean-tex:
	rm -f draco-autodoc.log draco-autodoc.aux
	rm -f draco-autodoc.toc draco-autodoc.lof draco-autodoc.lot

###############################################################################
#			 end of draco/doc/autodic Makefile
###############################################################################
