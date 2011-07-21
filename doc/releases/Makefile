###############################################################################
# draco/doc/releases Makefile
# Thomas M. Evans
# Thu Feb 11 10:32:47 1999
# Time-stamp: <04/04/28 12:37:08 tme>
###############################################################################
# @> Makefile for Draco build document
###############################################################################

DRACO := ../..
DRACOTEXINPUTS := ${DRACO}/environment/bibfiles
TEXINPUTSENVS := TEXINPUTS="${DRACOTEXINPUTS}:${TEXINPUTS}"

LATEXENVS := ${TEXINPUTSENVS}
LATEX := ${LATEXENVS} latex

BIBTEXENVS := BIBINPUTS="${DRACOTEXINPUTS}:${BSTINPUTS}"
BIBTEX := ${BIBTEXENVS} bibtex

MAKEIDX := makeindex

PDF := ps2pdf13

DVIPS := dvips
DVIPSFLAGS := -t letter -Ppdf

##---------------------------------------------------------------------------##
## SOURCES

tex_sources = draco-3_0_0.tex draco-3_2_0.tex draco-4_0_0.tex draco-5_0_0.tex
dvi_sources = $(patsubst %.tex,%.dvi,$(tex_sources))
ps_targets = $(patsubst %.dvi,%.ps,$(dvi_sources))
pdf_targets = $(patsubst %.dvi,%.pdf,$(dvi_sources))

##---------------------------------------------------------------------------##

.PHONY : all install clean pdf ps

all: $(dvi_sources)

install: 
	@echo install not ready

pdf: $(pdf_targets)

ps: $(ps_targets)

%.ps: %.dvi
	$(DVIPS) $(DVIPSFLAGS) $<

%.pdf: %.ps
	$(PDF) $<

%.dvi: %.tex
	$(LATEX) $(patsubst %.tex, %, $<)
	$(BIBTEX) $(patsubst %.tex, %, $<)
	$(LATEX) $(patsubst %.tex, %, $<)
	$(LATEX) $(patsubst %.tex, %, $<)

clean:
	-rm -f $(dvi_sources)
	-rm -f $(ps_targets)
	-rm -f $(pdf_targets)
	-rm -f *.blg
	-rm -f *.bbl
	-rm -f *.log
	-rm -f *.aux
	-rm -f *.toc *.lof *.lot
	-rm -f *~ *.flc

###############################################################################
#			 end of draco/doc/release Makefile
###############################################################################
