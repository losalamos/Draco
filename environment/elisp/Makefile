###############################################################################
#
# \file    draco/environment/elisp Makefile
# \author  Kelly G. Thompson
# \date    Mon May 23 10:33:22 2005
# \brief   Makefile for the Draco Elisp directory (compiles elisp
#          script to elc. files.)
#
###############################################################################
###############################################################################

elispsrc := $(wildcard *.el)
elcfiles := $(notdir $(patsubst %.el,%.elc,${elispsrc}))

.PHONY: all clean

all: ${elcfiles}

%.elc: %.el
	./compile-elisp $<

clean:
	rm -rf *.elc *~

##################################################
# Specific dependencies
##################################################

# doxymacs.elc: xml-parse.elc
# 	./compile-elisp $*.el

draco-config-modes.elc: draco-mode.elc draco-setup.elc draco-global-keys.elc draco-menu.elc
	./compile-elisp $*.el

draco-menu.elc: draco-mode.elc
	./compile-elisp $*.el

draco-new-files.elc: draco-helper-functions.elc
	./compile-elisp $*.el

draco-setup.elc: draco-faces.elc draco-new-files.elc
	./compile-elisp $*.el
