BIB_PATH = ..
CLS_PATH = ..

FILE = requirements

SOURCES = $(FILE).tex

%.pdf: %.ps
	ps2pdf -dEmbedAllFonts=true $<

$(FILE).ps: $(FILE).dvi
	dvips -P pdf -o $@ $<

$(FILE).dvi: $(SOURCES) 
	latex $(FILE)
	latex $(FILE)

level.dot: ../*.hh ../*.cc
	(cd ..;python ../../tools/include_tree.py)

level.ps: level.dot
	dot -Tps -o $@ $<

clean:
	-rm -rf $(FILE).dvi \
	        $(FILE).log *.aux texput.log \
		$(FILE).toc *.tex.flc \
                $(FILE).lof $(FILE).lot 
