hybridParticle.ps: hybridParticle.dvi
	dvips -o hybridParticle.ps hybridParticle.dvi

hybridParticle.dvi: hybridParticle.tex
	latex hybridParticle

clean:
	-rm hybridParticle.dvi hybridParticle.ps \
	        hybridParticle.log hybridParticle.aux 

