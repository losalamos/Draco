#  ========================================================================  #
# 
#       Author: Mark G. Gray
#               Los Alamos National Laboratory
#       Date:   Wed Apr  5 10:10:17 MDT 2000
# 
#       Copyright (c) 1999 U. S. Department of Energy. All rights reserved.
# 
#       $Id$
# 
#  ========================================================================  #

acf90.ps: acf90.dvi
	dvips -o acf90.ps acf90.dvi

acf90.dvi: acf90.tex
	latex acf90.tex
	bibtex acf90
	latex acf90.tex
	latex acf90.tex

view: acf90.ps
	ghostview acf90.ps

print: acf90.ps
	lpr acf90.ps

clean:
	rm -f *.aux *.log *.dvi *.ps *~

#  ========================================================================  #
