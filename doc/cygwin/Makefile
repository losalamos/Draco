###############################################################################
#
# \file    capsaicin/doc/vision_and_scope Makefile
# \author  Kelly G. Thompson
# \date    Tue Sept 19 10:33:22 2003
# \brief   Makefile for the Capsaicin Vision and Scope documents
#
###############################################################################
# $Id$
###############################################################################

latexsrc := $(wildcard *.tex)
psfiles := $(notdir $(patsubst %.tex,%.ps,${latexsrc}))
pdffiles := $(notdir $(patsubst %.tex,%.pdf,${latexsrc}))
dvifiles := $(notdir $(patsubst %.tex,%.dvi,${latexsrc}))

DRACO := ../..
DRACOTEXINPUTS := ${DRACO}/doc/tex
TEXINPUTSENVS := TEXINPUTS="${DRACOTEXINPUTS}:${TEXINPUTS}"

LATEXENVS := ${TEXINPUTSENVS}
LATEX := ${LATEXENVS} latex

BIBTEXENVS := BSTINPUTS="${DRACOTEXINPUTS}:${BSTINPUTS}"
BIBTEX := ${BIBTEXENVS} bibtex

MAKEIDX := makeindex

.PHONY : all xdvi ghostview install clean dvi postscript clean-tex index bib tex

all: postscript

#xdvi: dvi
#	xdvi ${FILE}

install: all
	@echo install not ready

postscript: ${psfiles}
	@echo "done creating postscript files: ${psfiles}"

dvi: ${dvifiles}

pdf: ${pdffiles}

%.ps: %.dvi
	dvips -o $@ $<

%.pdf: %.ps
	ps2pdf $<

%.dvi: %.tex 
	@echo ${LATEX} $<
	${LATEX} $<
#	@echo ${BIBTEX} $*
#	${BIBTEX} $*
#	@echo ${LATEX} $<
#	${LATEX} $<
#	@echo ${LATEX} $<
#	${LATEX} $<

index: 
	${MAKEIDX} ${latexsrc}

bib:
	${BIBTEX} ${latexsrc}

tex:
	${LATEX} ${latexsrc}

clean:
	rm -f *.dvi *.ps *.log *.aux *.pdf
	rm -f *.blg *.bbl
	rm -f *.idx *.ilg *.ind
	rm -f *.toc *.lof *.lot
	rm -f *.flc *~

clean-tex:
	rm -f *.log *.aux 
	rm -f *.toc *.lof *.lot

###############################################################################
# end of capsaicin/doc/vision_and_scope/Makefile
###############################################################################
