###############################################################################
# draco/doc/build Makefile
# Thomas M. Evans
# Thu Feb 11 10:32:47 1999
# Time-stamp: <03/05/09 17:01:01 tme>
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

.PHONY : all install clean postscript pdf real-clean

all: dbs.dvi 

install: 
	@echo install not ready

postscript: dbs.ps

pdf: dbs.pdf

%.dvi: %.tex $(source)
	${LATEX}  $(patsubst %.tex, %, $<)
	${BIBTEX} $(patsubst %.tex, %, $<)
	${MAKEIDX}$(patsubst %.tex, %, $<)
	${LATEX}  $(patsubst %.tex, %, $<)
	${LATEX}  $(patsubst %.tex, %, $<)

%.ps: %.dvi
	${DVIPS} ${DVIPSFLAGS} $<

%.pdf: %.ps
	${PDF} $<

clean:
	rm -f dbs.dvi dbs.log dbs.aux
	rm -f dbs.blg dbs.bbl
	rm -f dbs.idx dbs.ilg dbs.ind
	rm -f dbs.toc dbs.lof dbs.lot
	rm -f *~

real-clean: clean
	rm -f dbs.ps dbs.pdf

###############################################################################
#			 end of draco/doc/build Makefile
###############################################################################
