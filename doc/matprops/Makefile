MultiMaterialCell.ps: MultiMaterialCell.dvi
	dvips -o MultiMaterialCell.ps MultiMaterialCell.dvi

MultiMaterialCell.dvi: MultiMaterialCell.bbl MultiMaterialCell.tex
	latex MultiMaterialCell
	latex MultiMaterialCell

MultiMaterialCell.bbl: MultiMaterialCell.tex
	latex MultiMaterialCell
	bibtex MultiMaterialCell

clean:
	-rm MultiMaterialCell.dvi MultiMaterialCell.ps \
	        MultiMaterialCell.log MultiMaterialCell.aux \
		MultiMaterialCell.blg MultiMaterialCell.bbl
