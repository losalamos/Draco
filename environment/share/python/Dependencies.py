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
def get_files(dir):
    "Get a list of C/C++ files in the specified directory"
    
    import os,glob

    os.chdir(dir)

    files = glob.glob("*.hh")
    files.extend(glob.glob("*.cc"))
    files.extend(glob.glob("*.c"))
    files.extend(glob.glob("*.h"))

    return files

##---------------------------------------------------------------------------##
def file_includes(file):
    """
    Create a dictionary of dependencies by parsing 'file'. Keys in the
    dictionary are components and the values are filenames in that component.

       d:: component -> [list of included files from component]

    The included files appear in #include directives as either:

      #include \"package/filename\"   or
      #include <package/filename>

    For example:

    >>> file_includes('/home/mwbuksas/work/source/clubimc/head/src/imc/Source.hh')
    {'rng': ['Random.hh'], 'ds++': ['SP.hh'], 'mc': ['Particle_Stack.hh', 'Topology.hh']}
    
    """
    includes = {}
    f = open(file, 'r')

    # Loop through the file and look for include statements
    for line in f.readlines():

        # Check for "#include" from other packages
        match = quotes.match(line) or brackets.match(line)

        if match:

            package  = match.group('package')
            filename = match.group('filename')
            
            includes.setdefault(package, []).append(filename)

    f.close()

    return includes


##---------------------------------------------------------------------------##
def file_dependencies(file):
    """Extract just the components that contain header files that
    'file' depends on.

    >>> file_dependencies('/home/mwbuksas/work/source/clubimc/head/src/imc/Source.hh')
    ['rng', 'ds++', 'mc']

    """

    return file_includes(file).keys()



##---------------------------------------------------------------------------##
def directory_includes(directory):
    """Generates a combined depednency map for all of the files in a
    directory. 

    map :: component -> [list of included files]

    Examples:

    >>> d = directory_includes('/home/mwbuksas/work/source/clubimc/head/src/mc/')
    >>> print d['ds++']
    ['SP.hh', 'Assert.hh', 'Range_Finder.hh', 'Index_Converter.hh', 'Index_Counter.hh', 'Packing_Utils.hh', 'Soft_Equivalence.hh', 'Safe_Divide.hh']
    """

    import Utils

    dir_includes = {}

    for file in get_files(directory):
        for (component, files) in file_includes(file).items():
            Utils.unique_extend(dir_includes.setdefault(component, []), files)

    return dir_includes


##---------------------------------------------------------------------------##
def directory_dependencies(directory):
    """Extract the components that files in 'directory' depend on.

   >>> c = directory_dependencies('/home/mwbuksas/work/source/clubimc/head/src/imc/')
   >>> print c
   ['ds++', 'mc', 'cdi', 'rng', 'utils', 'c4', 'imc']
   """
    import Utils

    components = []
    for file in get_files(directory):
        Utils.unique_extend(components, file_dependencies(file))

    return components
    

##---------------------------------------------------------------------------##
def find_components_files(dependency_dict, exclude = []):

    """Extract a dictionary from 'dependency_dict' whose keys are
    referenced components and whose values are lists of included files.
    """

    import Utils

    components = {}

    for (file, include_dict) in dependency_dict.items():

        for (component, included_files) in include_dict.items():

            if component not in exclude:

                current_files = components.setdefault(component,[])
                Utils.unique_extend(current_files, included_files)

    return components

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


