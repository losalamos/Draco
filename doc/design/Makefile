###############################################################################
#
# \file    capsaicin/doc/design Makefile
# \author  Kelly G. Thompson
# \date    Tue Nov 16 8:20:05 2004
# \brief   Makefile for the Capsaicin dot levelization diagram
#
###############################################################################
# $Id$
###############################################################################

DOT=dot
DOT_OPTS=-Tps

PS2PDF=ps2pdf

dotsrc := $(wildcard *.dot)
psfiles := $(notdir $(patsubst %.dot,%.ps,${dotsrc}))
pdffiles := $(notdir $(patsubst %.dot,%.pdf,${dotsrc}))

all: ${psfiles}

pdf: ${pdffiles}

%.ps: %.dot
	${DOT} ${DOT_OPTS} -o $@ $<

%.pdf: %.ps
	${PS2PDF} $<

clean:
	rm -rf ${psfiles} ${pdffiles}
	rm -rf *~ *.flc

