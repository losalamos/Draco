###############################################################################
# draco/doc/build Makefile
# Thomas M. Evans
# Thu Feb 11 10:32:47 1999
# Time-stamp: <99/02/11 16:17:20 tme>
###############################################################################
# @> Makefile for Draco build document
###############################################################################

DRACO := ../..
DRACOTEXINPUTS := ${DRACO}/doc/tex
TEXINPUTSENVS := TEXINPUTS="${TEXINPUTS}:${DRACOTEXINPUTS}"

LATEXENVS := ${TEXINPUTSENVS}
LATEX := ${LATEXENVS} latex

BIBTEXENVS := BSTINPUTS="${BSTINPUTS}:${DRACOTEXINPUTS}"
BIBTEX := ${BIBTEXENVS} bibtex

MAKEIDX := makeindex

.PHONY : all install clean dvi postscript clean-tex index bib tex

all: postscript

install: all
	@echo install not ready

postscript: draco_bs.ps

dvi: draco_bs.dvi

draco_bs.ps: draco_bs.dvi
	dvips -o draco_bs.ps draco_bs.dvi

draco_bs.dvi: draco_bs.tex intro.tex model.tex compile.tex \
	      extern.tex adding.tex extend.tex draco_bs.sty
	${LATEX} draco_bs
	${BIBTEX} draco_bs
	${MAKEIDX} draco_bs
	${LATEX} draco_bs
	${LATEX} draco_bs

index: 
	${MAKEIDX} draco_bs

bib:
	${BIBTEX} draco_bs

tex:
	${LATEX} draco_bs

clean:
	rm -f draco_bs.dvi draco_bs.ps draco_bs.log draco_bs.aux
	rm -f draco_bs.blg draco_bs.bbl
	rm -f draco_bs.idx draco_bs.ilg draco_bs.ind
	rm -f draco_bs.toc draco_bs.lof draco_bs.lot

clean-tex:
	rm -f draco_bs.log draco_bs.aux
	rm -f draco_bs.toc draco_bs.lof draco_bs.lot

###############################################################################
#			 end of draco/doc/build Makefile
###############################################################################
