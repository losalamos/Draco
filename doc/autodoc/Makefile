draco-autodoc.ps: draco-autodoc.dvi
	dvips -o draco-autodoc.ps draco-autodoc.dvi

draco-autodoc.dvi: draco-autodoc.tex
	latex draco-autodoc
	bibtex draco-autodoc
	latex draco-autodoc
	latex draco-autodoc

clean:
	rm -f draco-autodoc.dvi draco-autodoc.ps \
	        draco-autodoc.log draco-autodoc.aux \
		draco-autodoc.bbl draco-autodoc.blg

