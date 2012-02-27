###############################################################################
# draco/doc/baltab Makefile
# John M. McGhee
# Time-stamp: <99/04/27 15:40:10 mcghee>
###############################################################################
# @> Makefile for Draco balance table document
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

postscript: draco-baltab.ps

dvi: draco-baltab.dvi

draco-baltab.ps: draco-baltab.dvi
	dvips -o draco-baltab.ps draco-baltab.dvi

draco-baltab.dvi: draco-baltab.tex 
	${LATEX} draco-baltab
#	${BIBTEX} draco-baltab
	${LATEX} draco-baltab
	${LATEX} draco-baltab

index: 
	${MAKEIDX} draco-baltab

bib:
	${BIBTEX} draco-baltab

tex:
	${LATEX} draco-baltab

clean:
	rm -f draco-baltab.dvi draco-baltab.ps draco-baltab.log draco-baltab.aux
	rm -f draco-baltab.blg draco-baltab.bbl
	rm -f draco-baltab.idx draco-baltab.ilg draco-baltab.ind
	rm -f draco-baltab.toc draco-baltab.lof draco-baltab.lot

clean-tex:
	rm -f draco-baltab.log draco-baltab.aux
	rm -f draco-baltab.toc draco-baltab.lof draco-baltab.lot

###############################################################################
#			 end of draco/doc/baltab Makefile
###############################################################################
