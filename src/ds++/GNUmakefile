ifdef nj
subnj = -j $(nj)
endif

default:
	$(MAKE) $(subnj) -f makefile incldeps=no deps
	$(MAKE) $(subnj) -f makefile

clean:
	$(MAKE) -f makefile clean
