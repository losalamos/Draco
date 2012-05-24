###############################################################################
#
# \file    draco/doc/gettingStarted Makefile
# \author  Kelly G. Thompson
# \date    Tue Apr 12 10:08:22 2000
# \brief   Makefile for the "Getting Started with Draco" technical memo.
#
###############################################################################
# $Id$
###############################################################################

FILE := gettingStarted

DRACO := ../..
DRACOTEXINPUTS := ${DRACO}/environment/latex
TEXINPUTSENVS := TEXINPUTS="${DRACOTEXINPUTS}:${TEXINPUTS}"

LATEX := ${TEXINPUTSENVS} latex
DVIPS := ${TEXINPUTSENVS} dvips

BIBTEXINPUTS := ${DRACO}/environment/bibtex
BIBTEXENVS := BSTINPUTS="${BIBTEXINPUTS}:${BSTINPUTS}"
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

pdf: postscript
	ps2pdf ${FILE}.ps

${FILE}.ps: ${FILE}.dvi
	${DVIPS} -o ${FILE}.ps ${FILE}.dvi

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
# end of draco/doc/gettingStarted/Makefile
###############################################################################
