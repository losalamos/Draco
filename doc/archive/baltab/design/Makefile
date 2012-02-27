###############################################################################
# draco/doc/baltab/design Makefile
# John M. McGhee
# Time-stamp: <99/10/01 15:39:23 mcghee>
###############################################################################
# @> Makefile for Draco balance table design document
###############################################################################

DRACO := ../../..
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

postscript: baltab-design.ps

dvi: baltab-design.dvi

baltab-design.ps: baltab-design.dvi
	dvips -o baltab-design.ps baltab-design.dvi

baltab-design.dvi: baltab-design.tex 
	${LATEX} baltab-design
	${BIBTEX} baltab-design
	${LATEX} baltab-design
	${LATEX} baltab-design

index: 
	${MAKEIDX} baltab-design

bib:
	${BIBTEX} baltab-design

tex:
	${LATEX} baltab-design

clean:
	rm -f baltab-design.dvi baltab-design.ps baltab-design.log baltab-design.aux
	rm -f baltab-design.blg baltab-design.bbl
	rm -f baltab-design.idx baltab-design.ilg baltab-design.ind
	rm -f baltab-design.toc baltab-design.lof baltab-design.lot

clean-tex:
	rm -f baltab-design.log baltab-design.aux
	rm -f baltab-design.toc baltab-design.lof baltab-design.lot

###############################################################################
#			 end of draco/doc/baltab/design Makefile
###############################################################################
