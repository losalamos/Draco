###############################################################################
# draco/doc/build Makefile
# Makefile for Draco build document
# $Id$
###############################################################################

DRACO := ../..
DRACOTEXINPUTS := ${DRACO}/environment/latex
TEXINPUTSENVS := TEXINPUTS="${DRACOTEXINPUTS}:${TEXINPUTS}"

LATEXENVS := ${TEXINPUTSENVS}
LATEX := ${LATEXENVS} latex

BIBTEXENVS := BSTINPUTS="${DRACOTEXINPUTS}:${BSTINPUTS}"
BIBTEX := ${BIBTEXENVS} bibtex

MAKEIDX := makeindex

DVIPS := dvips
DVIPSFLAGS := -t letter -Ppdf

PDF := ps2pdf13

##---------------------------------------------------------------------------##
## DEPENDENCIES
##---------------------------------------------------------------------------##

texsource := $(wildcard *.tex)
figsource := $(wildcard *.eps)
source    := $(texsource) $(figsource)

##---------------------------------------------------------------------------##
## RULES
##---------------------------------------------------------------------------##

.PHONY : all install clean dvi postscript index bib tex pdf 

all: dbs.pdf

install: 
	@echo install not ready

postscript: dbs.ps

pdf: dbs.pdf

%.ps: %.dvi
	${DVIPS} ${DVIPSFLAGS} $<

%.pdf: %.ps
	${PDF} $<

%.dvi: %.tex $(source)
	${LATEX}  $(patsubst %.tex, %, $<)
	${BIBTEX} $(patsubst %.tex, %, $<)
	${MAKEIDX}$(patsubst %.tex, %, $<)
	${LATEX}  $(patsubst %.tex, %, $<)
	${LATEX}  $(patsubst %.tex, %, $<)

index: 
	${MAKEIDX} dbs

bib:
	${BIBTEX} dbs

tex:
	${LATEX} dbs

clean:
	rm -f dbs.dvi dbs.ps dbs.log dbs.aux
	rm -f dbs.blg dbs.bbl dbs.pdf
	rm -f dbs.idx dbs.ilg dbs.ind
	rm -f dbs.toc dbs.lof dbs.lot
	rm -f *~ dbs.ps

###############################################################################
#			 end of draco/doc/build Makefile
###############################################################################
