radTransfer.ps: radTransfer.dvi
	dvips -o radTransfer.ps radTransfer.dvi

radTransfer.dvi: radTransfer.bbl radTransfer.tex
	latex radTransfer
	latex radTransfer

radTransfer.bbl: radTransfer.tex
	latex radTransfer
	bibtex radTransfer

clean:
	-rm radTransfer.dvi radTransfer.ps radTransfer.log radTransfer.aux \
		radTransfer.blg radTransfer.bbl
