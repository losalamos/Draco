#!/usr/bin/env python
#======================================================================
# module: Packages
#
# Contains functions for parsing C/C++ code, directories and packages
# for tasks such as component dependency tracking.
#
# Michael W. Buksas
#
# $Id$
#======================================================================


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
def get_dependencies(file):

    import re

    """

    Create a dictionary of dependencies by parsing
    'file'. Viewed as a functions, the resulting signature is:

       get_dependencies:: file -> d, where:

       d:: package -> [list of included files]

    The included files appear in #include directives as either:

      #include \"package/filename\"   or
      #include <package/filename>

    >>> get_dependencies('/codes/radtran/vendors/clubimc/clubimc/src/imc/Source.hh')
    {'rng': ['Random.hh'], 'ds++': ['SP.hh'], 'mc': ['Particle_Stack.hh', 'Topology.hh']}
    
    """

    include_dict = {}

    include  = "#include\s*"
    contents = "(?P<package>[\w+]*)\/+(?P<filename>[\w+.]*.\w*)\s*"

    quotes   = re.compile(include + '\"'+ contents + '\"')
    brackets = re.compile(include + '<' + contents + '>')
    
    f = open(file, 'r')

    # Loop through the file and get look for include statements
    for line in f.readlines():

        # Check for "#include" from other packages
        match = quotes.match(line) or brackets.match(line)

        if match:

            package  = match.group('package')
            filename = match.group('filename')
            
            include_dict.setdefault(package, []).append(filename)

    f.close()

    return include_dict

##---------------------------------------------------------------------------##
def scan_directory(directory):

    """Build a dictionary which maps filenames in the given directory 
    to the output of get_dependencies(filename). Viewing the resulting
    dictionary as a map, it's signature would be:

       d: filename -> package -> [included_files]

   >>> d = scan_directory('/codes/radtran/vendors/clubimc/clubimc/src/imc/')

   >>> d['Source.hh']
   {'rng': ['Random.hh'], 'ds++': ['SP.hh'], 'mc': ['Particle_Stack.hh', 'Topology.hh']}

   >>> d['Global.cc']
   {}

   >>> d['Global.hh']
   {'mc': ['Math.hh', 'Constants.hh']}
    
    """
    
    files = get_files(directory)
    file_includes = {}

    for file in files:
        file_includes[file] = get_dependencies(file)

    return file_includes



##---------------------------------------------------------------------------##
def find_components(dependency_dict, exclude = []):

    """ Extract a list of referenced components from 'dependency_dict',
    excluding the ones in list 'exclude'.

   >>> d = dependency_dict('/codes/radtran/vendors/clubimc/clubimc/src/imc/')
    
   >>> c = find_components(d)
   >>> print c
   ['imc', 'mc', 'c4', 'rng', 'cdi', 'ds++']

   >>> d = dependency_dict('/codes/radtran/vendors/clubimc/clubimc/src/imc/test/')
   >>> find_components(d, c)
   ['cdi_analytic']
    
    """

    import Utils

    components = []

    for (file, include_dict) in dependency_dict.items():

        for component in include_dict.keys():

            if (component not in exclude):
                Utils.unique_append(components, component)

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
    import doctest, Pacakges
    doctest.testmod(Packages)


##---------------------------------------------------------------------------##
## Main Program
##---------------------------------------------------------------------------##

if __name__=='__main__':

    """Module Packages

    Contains functions useful for pasrsing the contents of C/C++
    code, and determining dependencies between files and components.

    See the help information for individual functions for more
    information. 

    """

    # Run the unit tests
    _test()


