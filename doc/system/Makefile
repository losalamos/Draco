###############################################################################
#
# \file    draco/doc/system/Makefile
# \brief   Makefile for rn98046.tex
#
###############################################################################
# $Id$
###############################################################################

FILE := rn98046

DRACO := ../..
DRACOTEXINPUTS := ${DRACO}/doc/tex
TEXINPUTSENVS := TEXINPUTS="${DRACOTEXINPUTS}:${TEXINPUTS}"

LATEXENVS := ${TEXINPUTSENVS}
LATEX := ${LATEXENVS} latex

DRACOBIBDIR=${DRACO}/environment/bibfiles
BIBTEXENVS := BSTINPUTS="${DRACOBIBDIR}:${BSTINPUTS}"

BIBTEX := ${BIBTEXENVS} bibtex 
MAKEIDX := makeindex

.PHONY : all xdvi ghostview install clean dvi postscript clean-tex index bib tex

all: postscript

ghostview: postscript
	ghostview ${FILE}.ps

xdvi: dvi
	xdvi ${FILE}

install: all
	@echo install not ready

postscript: ${FILE}.ps

dvi: ${FILE}.dvi

pdf: ${FILE}.pdf

${FILE}.pdf: ${FILE}.ps
	ps2pdf ${FILE}.ps

${FILE}.ps: ${FILE}.dvi
	dvips -o ${FILE}.ps ${FILE}.dvi

${FILE}.dvi: ${FILE}.tex 
	${LATEX} ${FILE}
	${BIBTEX} ${FILE}
	${LATEX} ${FILE}
	${LATEX} ${FILE}

index: 
	${MAKEIDX} ${FILE}

bib:
	${BIBTEX} ${FILE}

tex:
	${LATEX} ${FILE}

clean:
	rm -f ${FILE}.dvi ${FILE}.ps ${FILE}.log ${FILE}.aux
	rm -f ${FILE}.blg ${FILE}.bbl
	rm -f ${FILE}.idx ${FILE}.ilg ${FILE}.ind
	rm -f ${FILE}.toc ${FILE}.lof ${FILE}.lot
	rm -f *.flc *~

clean-tex:
	rm -f ${FILE}.log ${FILE}.aux
	rm -f ${FILE}.toc ${FILE}.lof ${FILE}.lot

###############################################################################
# end of draco/doc/system/Makefile
###############################################################################
