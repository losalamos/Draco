###############################################################################
# python script useful for updating milagro's benchmark output files for
# regression.  the script merely copies *.test to *.bench for all the 
# regression tests.  the user can modify the paths to copy the *.test in
# build directory directly to the *.bench in the source directory, if desired.
#
# Todd Urbatsch, CCS-4, LANL, 25 SEP 2001
#
###############################################################################

import os, string

parallel_descriptor = [("inf",34),     ("str",31),    ("tp",7), \
                       ("p2_inf",4),   ("p2_str",6),  ("start",6), \
                       ("restart",12), ("b_start",1), ("b_restart",2)]
serial_descriptor   = [("inf",34),     ("str",31),    ("tp",7), \
                       ("start",3),    ("restart",3), \
                       ("b_start",1),  ("b_restart",1)]

# number: the two digit character strings necessary to build filenames.
number = ["01","02","03","04","05","06","07","08","09","10",
          "11","12","13","14","15","16","17","18","19","20",
          "21","22","23","24","25","26","27","28","29","30",
          "31","32","33","34","35","36","37","38","39","40",
          "41","42","43","44","45","46","47","48","49","50"]

descriptor = serial_descriptor

group_num = 0
for desc in descriptor:
    group     = desc[0]
    num_cases = desc[1]
    case = 0
    while case < num_cases:
        name  = group + number[case]
        cpfrom = name + ".test"
        cpto   = name + ".bench"
        print "mv -f %s %s" % (cpfrom, cpto)
        os.system("cp %s %s" % (cpfrom, cpto))
        case = case + 1
    group_num = group_num + 1

