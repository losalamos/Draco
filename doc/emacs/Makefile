###############################################################################
# draco/doc/emacs Makefile
# Thomas M. Evans
# Thu Feb 11 10:32:47 1999
# Time-stamp: <99/02/22 11:21:24 tme>
###############################################################################
# @> Makefile for Draco build document
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

postscript: draco-emacs.ps

dvi: draco-emacs.dvi

draco-emacs.ps: draco-emacs.dvi
	dvips -o draco-emacs.ps draco-emacs.dvi

draco-emacs.dvi: draco-emacs.tex 
	${LATEX} draco-emacs
	${BIBTEX} draco-emacs
	${LATEX} draco-emacs
	${LATEX} draco-emacs

index: 
	${MAKEIDX} draco-emacs

bib:
	${BIBTEX} draco-emacs

tex:
	${LATEX} draco-emacs

clean:
	rm -f draco-emacs.dvi draco-emacs.ps draco-emacs.log draco-emacs.aux
	rm -f draco-emacs.blg draco-emacs.bbl
	rm -f draco-emacs.idx draco-emacs.ilg draco-emacs.ind
	rm -f draco-emacs.toc draco-emacs.lof draco-emacs.lot

clean-tex:
	rm -f draco-emacs.log draco-emacs.aux
	rm -f draco-emacs.toc draco-emacs.lof draco-emacs.lot

###############################################################################
#			 end of draco/doc/emacs Makefile
###############################################################################
