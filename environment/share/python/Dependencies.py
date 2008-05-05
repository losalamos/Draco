#!/usr/bin/env python
#======================================================================
# module: Dependencies
#
# Contains functions for parsing C/C++ code, directories and packages
# for tasks such as component dependency tracking.
#
# Michael W. Buksas
#
# $Id$
#======================================================================

import re

# Variables for #include regular expressions.
include  = "#include\s*"
contents = "(?P<package>[\w+]*)\/+(?P<filename>[\w+.]*.\w*)\s*"

quotes   = re.compile(include + '\"'+ contents + '\"')
brackets = re.compile(include + '<' + contents + '>')
    

##---------------------------------------------------------------------------##
def gen_cpp_files(dir, recurse=True):
    """Generate a list of C/C++ files in the specified
    directory. Recurse into subdirectories if recurse=True."""

    import os
    import fnmatch

    extensions = ["hh", "cc", "h", "c"]

    for path, dirlist, filelist in os.walk(dir):
        if not recurse: dirlist[:] = []
        for ext in extensions:
            for name in fnmatch.filter(filelist, "*.%s" % ext):
                yield os.path.join(path,name)


##---------------------------------------------------------------------------##
def get_files(dir, recurse=True):
    import FileUtils
    return FileUtils.gen_open(gen_cpp_files(dir, recurse))

##---------------------------------------------------------------------------##
def file_includes(files):
    """
    Create a dictionary of dependencies by parsing each file in
    'files'. Keys in the dictionary are components imported from and
    the values are the names of imported files.

    file -> {component: set of files included into file}

    The included files appear in #include directives as either:

      #include \"package/filename\"   or
      #include <package/filename>

    For example:

    >>> file_includes([open('/home/mwbuksas/work/source/clubimc/head/src/imc/Source.hh','r')])
    {'rng': set(['Random.hh']), 'ds++': set(['SP.hh']), 'mc': set(['Topology.hh', 'Particle_Stack.hh'])}
    
    """

    import itertools

    includes = {}
    for line in itertools.chain(*files):

        match = quotes.match(line) or brackets.match(line)

        if match:
            package  = match.group('package')
            filename = match.group('filename')
            includes.setdefault(package, set()).add(filename)

    return includes




##---------------------------------------------------------------------------##
## Reducers:
##---------------------------------------------------------------------------##
# These pair down the information from file_includes:

##---------------------------------------------------------------------------##
def file_inc_comps(files):
    """Extract just the components that contain header files that
    'file' depends on.

    >>> file_inc_comps([open('/home/mwbuksas/work/source/clubimc/head/src/imc/Source.hh','r')])
    set(['rng', 'ds++', 'mc'])

    """

    return set(file_includes(files).keys())


##---------------------------------------------------------------------------##
def dir_inc_file(directory):
    """Generates a merged depednency map for all of the files in a
    directory. 

    map :: directory -> {set of included files}

    >>> d = dir_inc_file('/home/mwbuksas/work/source/clubimc/head/src/mc/')
    >>> print d['ds++']
    set(['Index_Counter.hh', 'Safe_Divide.hh', 'Assert.hh', 'SP.hh', 'Soft_Equivalence.hh', 'Range_Finder.hh', 'Packing_Utils.hh', 'Index_Converter.hh'])
    """

    return file_includes(get_files(directory, False))

##---------------------------------------------------------------------------##
def dir_inc_comp(directory):
    """Extract the components that files in 'directory' include from.

   >>> dir_inc_comp('/home/mwbuksas/work/source/clubimc/head/src/imc/')
   set(['imc', 'mc', 'utils', 'rng', 'cdi', 'ds++', 'c4'])
   """
    return file_inc_comps(get_files(directory, False))



##---------------------------------------------------------------------------##
## Test function
##---------------------------------------------------------------------------##
def _test():
    import doctest, Dependencies
    doctest.testmod(Dependencies)


##---------------------------------------------------------------------------##
## Main Program
##---------------------------------------------------------------------------##

if __name__=='__main__':

    """Module Dependencies

    Contains functions useful for pasrsing the contents of C/C++
    code, and determining dependencies between files and components.

    See the help information for individual functions for more
    information. 

    """

    # Run the unit tests
    _test()


