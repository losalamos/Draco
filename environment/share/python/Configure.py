"""A module to support operations related to running configure
scripts.

"""


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


