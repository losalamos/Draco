###############################################################################
# draco/doc/GenericProgInSolon Makefile
# Randy M. Roberts
# Time-stamp: <99/02/11 08:53:01 rsqrd>
###############################################################################
# @> Makefile for Generic Programming in the SOLON Interface (U)
###############################################################################

DRACO := ../..
DRACOTEXINPUTS := ${DRACO}/doc/tex
TEXINPUTSENVS := TEXINPUTS="${TEXINPUTS}:${DRACOTEXINPUTS}"

LATEXENVS := ${TEXINPUTSENVS}
LATEX := ${LATEXENVS} latex

BIBTEXENVS := ${TEXNPUTSENVS}
BIBTEX := ${BIBTEXENVS} bibtex

.PHONY : all install clean dvi postscript

all: postscript

install: all
	@echo install not ready

postscript: Solon.ps

dvi: Solon.dvi

Solon.ps: Solon.dvi
	dvips -o Solon.ps Solon.dvi

Solon.dvi: Solon.tex Solonpreamble.tex mesh.uml.eps p13t.uml.eps \
   diff.uml.eps Solon.bib
	${LATEX} Solon
	${BIBTEX} Solon
	${LATEX} Solon
	${LATEX} Solon

clean:
	-rm Solon.dvi Solon.ps Solon.log Solon.aux Solon.bbl Solon.blg
