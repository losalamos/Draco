"""A module to support operations related to running configure
scripts.

"""

kind_options = {"opt"   : "--with-dbc=0 --with-opt=3 --disable-debug",
                "debug" : "--with-dbc=7 --enable-debug"
                }

known_kinds = kind_options.keys()

compiler_names = {"pgCC" : "pgi",
                  "g++"  : "gcc"
                  }

##---------------------------------------------------------------------------##
def make_option_string(name, lib='', inc=''):
    """Make a configure string from a dictionary object containg it's name
    and path locations. The resulting string is of the form:
    
    --with-<name>-lib=<lib> --with-<name>-inc=<inc>
    
    where the unspecified parts come from the dictionary argument values
    with those keys. Either part of the string will be omitted if the
    corresponding key is not present in the dictionary.
    """

    parts = []
    if lib:
        parts.append("--with-%s-lib=%s" % (name, lib))
    if inc:
        parts.append("--with-%s-inc=%s" % (name, inc))

    return " ".join(parts)


##---------------------------------------------------------------------------##
def make_vendor_string(vendor_libs):
    """Make a configure string for the provided list of vendors."""

    vendor_string = " ".join(
        [Configure.make_option_string(**platform.vendors[vendor]) 
         for vendor in vendor_libs])

    return vendor_string
    
##---------------------------------------------------------------------------##
def get_compiler_string():
    """Determine the desired compiler from the environment variable CXX, if
    possible, and return a corresponding configuration string. If the
    environment variable isn't set, return the empty string."""

    environ_cc = os.path.basename(os.environ.get("CXX", ""))

    compiler = compiler_names.get(environ_cc)
    if compiler: return "--with-cxx=%s" % compiler
    else: return ""


##---------------------------------------------------------------------------##
def get_kinds(given_kinds):
    """Expand the given kinds names into complete kinds names. Throws
    an exception if any names can't be resolved."""

    # Expand the "kinds" arguments.
    try:
        kinds = [Utils.disambiguate(k, known_kinds) for k in
                 options.build_kinds]

    except KeyError, e:
        die("Could not understand kind argument: %s" % e[0])

    return kinds



##---------------------------------------------------------------------------##
def get_mpi_string(options):

    # MPI: Use the module system to provide values.
    if options.serial:
        return ""
    else:
        mpi_inc_dir = os.environ.get("MPI_INC_DIR") or die("Variable MPI_INC_DIR not set")
        mpi_lib_dir = os.environ.get("MPI_LIB_DIR") or die("Variable MPI_LIB_DIR not set")

    return " ".join([
        "--with-c4=mpi",
        "--with-mpi=lampi",
        "--with-mpi-inc=%s" % mpi_inc_dir,
        "--with-mpi-lib=%s" % mpi_lib_dir])



