###############################################################################
# draco/src/makefile
# Geoffrey Furnish
# 5 June 1997
###############################################################################
# @> Rules for building Draco components.
###############################################################################

export dodeps := no

ifndef nj
nj	:= 1
endif

subnj	:= -j $(nj)

all:
	(cd ds++;    $(MAKE) $(subnj))
	(cd nml;     $(MAKE) $(subnj))
	(cd c4;      $(MAKE) $(subnj))
	(cd geom;    $(MAKE) $(subnj))
	(cd mesh;    $(MAKE) $(subnj))
	(cd draco;   $(MAKE) $(subnj))

t:
	(cd ds++;    $(MAKE) $(subnj) incldeps=no deps; $(MAKE) $(subnj))
	(cd nml;     $(MAKE) $(subnj) incldeps=no deps; $(MAKE) $(subnj))
	(cd c4;      $(MAKE) $(subnj) incldeps=no deps; $(MAKE) $(subnj))

some:
	(cd ds++;    $(MAKE) $(subnj))
	(cd nml;     $(MAKE) $(subnj))
	(cd c4;      $(MAKE) $(subnj))

notests:
	(cd ds++;    $(MAKE))
	(cd nml;     $(MAKE) testprogs=no)
	(cd c4;      $(MAKE))

cleanall: clean
	$(RM) *~
	$(RM) -rf *repository
	$(RM) */*~
	$(RM) */.nml* */*_prs.cc

clean:
	(cd ds++;    $(MAKE) clean cleanonly=yes)
	(cd nml;     $(MAKE) clean cleanonly=yes)
	(cd c4;      $(MAKE) clean cleanonly=yes)
	(cd util;    $(MAKE) clean cleanonly=yes)
	(cd geom;    $(MAKE) clean cleanonly=yes)
	(cd mesh;    $(MAKE) clean cleanonly=yes)
	(cd draco;   $(MAKE) clean cleanonly=yes)

###############################################################################
#			   end of draco/src/makefile
###############################################################################
