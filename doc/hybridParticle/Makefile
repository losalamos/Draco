hybridParticle.ps: hybridParticle.dvi plot.eps
	dvips -o hybridParticle.ps hybridParticle.dvi

hybridParticle.dvi: hybridParticle.tex
	latex hybridParticle
	latex hybridParticle

clean:
	-rm hybridParticle.dvi hybridParticle.ps \
	        hybridParticle.log hybridParticle.aux 

