###############################################################################
# draco/doc/emacs Makefile
# Thomas M. Evans
# Thu Feb 11 10:32:47 1999
# Time-stamp: <99/02/22 11:21:24 tme>
###############################################################################
# @> Makefile for Draco build document
###############################################################################

FILE := draco-emacs

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

ghostview: postscript
	ghostview ${FILE}.ps

xdvi: dvi
	xdvi ${FILE}

pdf: postscript
	ps2pdf ${FILE}.ps

install: all
	@echo install not ready

postscript: ${FILE}.ps

dvi: ${FILE}.dvi

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
	@rm -f ${FILE}.dvi ${FILE}.ps ${FILE}.log ${FILE}.aux
	@rm -f ${FILE}.blg ${FILE}.bbl ${FILE}.pdf
	@rm -f ${FILE}.idx ${FILE}.ilg ${FILE}.ind
	@rm -f ${FILE}.toc ${FILE}.lof ${FILE}.lot
	@rm -f *~

clean-tex:
	rm -f ${FILE}.log ${FILE}.aux
	rm -f ${FILE}.toc ${FILE}.lof ${FILE}.lot

###############################################################################
#			 end of draco/doc/emacs Makefile
###############################################################################
