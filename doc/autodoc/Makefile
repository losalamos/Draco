draco-autodoc.ps: draco-autodoc.dvi
	dvips -o draco-autodoc.ps draco-autodoc.dvi

draco-autodoc.dvi: draco-autodoc.tex
	latex draco-autodoc

clean:
	-rm draco-autodoc.dvi draco-autodoc.ps \
	        draco-autodoc.log draco-autodoc.aux 

