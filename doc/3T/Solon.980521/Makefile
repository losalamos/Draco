solon.ps: solon.dvi
	dvips -o solon.ps solon.dvi

solon.dvi: SphereDoom.ps solon.tex solonpreamble.tex
	latex solon
	latex solon

clean:
	-rm solon.dvi solon.ps solon.log solon.aux
