###############################################################################
# /draco/doc/imctest/makefile 
# Thomas M. Evans
# Mon Jun  8 10:50:54 1998
###############################################################################
# @> Makefile for IMCTEST development Manual
###############################################################################

## MACROS

rm = rm -rf

## TeX srcs

texsrc = $(wildcard *.tex)
bibsrc = $(wildcard *.bib)
figsrc = $(wildcard *.cfg *.eps)

imc-dev.ps: imc-dev.dvi
	dvips imc-dev

imc-dev.dvi: $(texsrc) $(figsrc) imc-dev.bbl
	latex imc-dev
	latex imc-dev

imc-dev.bbl: $(texsrc) $(bibsrc)
	latex imc-dev
	bibtex imc-dev

clean: 
	$(rm) imc-dev.ps imc-dev.bbl imc-dev.dvi
	$(rm) imc-dev.aux imc-dev.blg imc-dev.log

###############################################################################
#			 end of /draco/doc/imctest/makefile
###############################################################################
