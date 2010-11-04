###############################################################################
# draco/doc/release Makefile
# Thomas M. Evans
# Thu Feb 11 10:32:47 1999
# Time-stamp: <99/05/05 17:00:33 tme>
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

.PHONY : all install clean dvi postscript clean-tex index bib tex pdf

all: postscript
pdf: draco-release.pdf

install: all
	@echo install not ready

postscript: draco-release.ps

dvi: draco-release.dvi

draco-release.ps: draco-release.dvi
	dvips -t letter -Ppdf -o draco-release.ps draco-release.dvi

draco-release.pdf: draco-release.ps
	ps2pdf $<

draco-release.dvi: draco-release.tex 
	${LATEX} draco-release
	${BIBTEX} draco-release
	${LATEX} draco-release
	${LATEX} draco-release

index: 
	${MAKEIDX} draco-release

bib:
	${BIBTEX} draco-release

tex:
	${LATEX} draco-release

clean:
	rm -f draco-release.dvi draco-release.ps draco-release.log draco-release.aux
	rm -f draco-release.blg draco-release.bbl
	rm -f draco-release.idx draco-release.ilg draco-release.ind
	rm -f draco-release.toc draco-release.lof draco-release.lot

clean-tex:
	rm -f draco-release.log draco-release.aux
	rm -f draco-release.toc draco-release.lof draco-release.lot

###############################################################################
#			 end of draco/doc/release Makefile
###############################################################################
