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
	(cd spdf;    $(MAKE) $(subnj))
	(cd util;    $(MAKE) $(subnj))
	(cd linalg;  $(MAKE) $(subnj))
	(cd geom;    $(MAKE) $(subnj))
	(cd mesh;    $(MAKE) $(subnj))
	(cd 3T;      $(MAKE) $(subnj))
	(cd sn;      $(MAKE) $(subnj))
	(cd imctest; $(MAKE) $(subnj))
#	(cd draco;   $(MAKE) $(subnj))

clean:
	(cd ds++;    $(MAKE) clean cleanonly=yes)
	(cd nml;     $(MAKE) clean cleanonly=yes)
	(cd c4;      $(MAKE) clean cleanonly=yes)
	(cd spdf;    $(MAKE) clean cleanonly=yes)
	(cd util;    $(MAKE) clean cleanonly=yes)
	(cd linalg;  $(MAKE) clean cleanonly=yes)
	(cd geom;    $(MAKE) clean cleanonly=yes)
	(cd mesh;    $(MAKE) clean cleanonly=yes)
	(cd 3T;      $(MAKE) clean cleanonly=yes)
	(cd sn;      $(MAKE) clean cleanonly=yes)
	(cd imctest; $(MAKE) clean cleanonly=yes)
	(cd pycxx;   $(MAKE) clean cleanonly=yes)
	(cd pydraco; $(MAKE) clean cleanonly=yes)
#	(cd draco;   $(MAKE) clean cleanonly=yes)

# The clean target has a pernicious failure mode.  If one of these
# directories does not exist (like you haven't checked it out), then
# the 'cd' will fail, leaving you running "make clean" in the current
# directory, so that the process loops an infinite number of times.
# Need to fix this!

###############################################################################
#			   end of draco/src/makefile
###############################################################################
