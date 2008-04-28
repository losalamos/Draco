"""A module for setting up the configuration and building of Jayenne
codes.
""" 

import Repo, CVS

class PackageError(Exception):
    """Reports errors relating to the Jayenne packages"""
    pass

class LookupError(PackageError): pass


# A map from package names to repositiry names.
repositories = {'draco'       : "draco",
                'tools'       : "draco",
                'imcdoc'      : "jayenne",
                'clubimc'     : "jayenne",
                'milagro'     : "jayenne",
                'wedgehog'    : "jayenne",
                'radtest'     : "jayenne",
                'jayenne'     : "jayenne",
                'uncleMcFlux' : "jayenne"}

packages     = repositories.keys()


##---------------------------------------------------------------------------##
## Respository and CVS source packages:
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
def is_valid_package(package): return package in packages
def disambiguate_package(package): return disambiguate(package, packages)

##---------------------------------------------------------------------------##
def get_repository(package_name):
    """Get a Repo.Repository object for the repository that contains the
    given package.

    >>> r = get_repository('uncleMcFlux')
    >>> r.name
    'jayenne'

    >>> r = get_repository('clubimc/src/mc')
    >>> r.name
    'jayenne'

    >>> r = get_repository('draco/environment')
    >>> r.name
    'draco'

    """
    parent = package_name.split("/",1)[0]
    assert(is_valid_package(parent))
    return Repo.get_repository(repositories[parent])


##---------------------------------------------------------------------------##
def make_package(package_name):
    """Get a CVS.Package object representing a Jayenne package.

    >>> w = make_package('wedgehog')
    >>> w.name
    'wedgehog'
    >>> w.repository.name
    'jayenne'
    >>> w.repository.location
    '/ccs/codes/radtran/cvsroot'
    >>> w.path
    '/ccs/codes/radtran/cvsroot/wedgehog'

    """

    repository = get_repository(package_name)
    return CVS.Package(package_name, repository)


##---------------------------------------------------------------------------##
def make_working_copy(package_name, checkout_name, tag_kind, tag_name):
    """Create a CVS.working_copy_object from a package name, the name
    to check it out under, and the tag information
    """

    tag = CVS.make_tag(tag_kind, tag_name)
    package = make_package(package_name)

    return CVS.WorkingCopy(package, tag, checkout_name)




##---------------------------------------------------------------------------##
def get_draco(all_deps, head = False, import_only = False):
    """A simple interface to the get_draco scripts which live in the
    top directory of Jayenne pacakges

    >>> get_draco(["draco", "clubimc"])
    './get_draco -Q -l -d /ccs/codes/radtran/cvsroot -j /ccs/codes/radtran/cvsroot'
    """

    command = "./get_draco -Q"

    if head: command += " -H"
    else:    command += " -l"

    if import_only: command += " -i"

    all_repos     = [get_repository(dep) for dep in all_deps]
    repos_by_name = [(repo.name, repo) for repo in all_repos]
    for name, repo in repos_by_name:
        if name == "draco":   command += " -d %s" % repo.location
        if name == "jayenne": command += " -j %s" % repo.location

    return command


##---------------------------------------------------------------------------##
def autoconf(name): 
    """The script which runs autoconf in Jayenne packages
    >>> autoconf('FooBar')
    ./FooBar_config

    """
    
    return "./%s_config" % name


##---------------------------------------------------------------------------##
## Dependencies and Components
##---------------------------------------------------------------------------##

# These are package level dependencies.
dependencies = {'draco'    : [],
                'clubimc'  : ['draco'],
                'milagro'  : ['draco', 'clubimc'],
                'wedgehog' : ['draco', 'clubimc']
                }

components = {'draco'    : ['RTT_Format_Reader', 
                            'c4', 
                            'cdi',
                            'cdi_analytic', 
                            'ds++', 
                            'meshReaders',
                            'mesh_element', 
                            'rng', 
                            'traits', 
                            'viz'
                           ],
              'clubimc'  : ['mc', 
                            'imc', 
                            'chimpy', 
                            'rng_nr', 
                            'utils'
                           ],
              'milagro'  : ['milagro',
                            'milagro_amr_rz',
                            'milagro_amr_rz_rz_mg',
                            'milagro_amr_xyz',
                            'milagro_amr_xyz_mg',
                            'milagro_builders',
                            'milagro_data',
                            'milagro_interfaces',
                            'milagro_manager',
                            'milagro_r',
                            'milagro_r_mg',
                            'milagro_release',
                            'milagro_rz',
                            'milagro_rz_mg',
                            'milagro_xyz',
                            'milagro_xyz_mg'
                            ],
              'wedgehog' : ['wedgehog',
                            'wedgehog_components',
                            'wedgehog_dd',
                            'wedgehog_gs',
                            'wedgehog_interfaces',
                            'wedgehog_managers',
                            'wedgehog_output',
                            'wedgehog_release',
                            'wedgehog_shunt',
                            'fortran_shunts'
                            ]
              }


##---------------------------------------------------------------------------##
def component_to_package(component):
    """Return the package that a component belongs to

    >>> component_to_package('c4')
    'draco'

    >>> component_to_package('bite me')
    Traceback (most recent call last):
    ...
    LookupError: Found no packages for bite me
    """
    results = [package for (package,comp_list) in components.items() if
               component in comp_list]

    if len(results) > 1: 
        raise LookupError("Found multiple packages for %s" % component)

    if len(results) < 1: 
        raise LookupError("Found no packages for %s" % component)
    
    return results[0]

##---------------------------------------------------------------------------##
def get_dependencies(package):
    """Get a list of packages the given package depends on.

    >>> get_dependencies('wedgehog')
    ['draco', 'clubimc']

    """
    assert(is_valid_package(package))
    return dependencies[package]


##---------------------------------------------------------------------------##
def is_install_dir(path):
    """Check to see if the given path looks like a install location for
    a Jayenne component"""

    return os.path.exists(os.path.join(path, "lib")) and os.path.exists(
        os.path.join(path, "include"))


##---------------------------------------------------------------------------##
def extract_jayenne_dependencies(target, options):

    """Return a map of the form: {component : path}

    Keys and values are both strings. The component keys are the
    jayenne components that 'target depends on. The path values may
    contain keywords for later textual substitution: <install>, <name>
    and <kind>.

    """

    assert(target in packages)

    depends = dependencies[target]

    depend_mapping = {}
    for component in depends:

        directory = getattr(options, component)
        depend_mapping[component] = directory

    return depend_mapping
    

##---------------------------------------------------------------------------##
def expand_template(path, keywords):
    """Expand a string with <foobar> keywords in it using the values
    given in the keyword arguments

    >>> expand_template("<install>/really/<kind>",\
                         dict(install="/some/place", kind="nice"))
    '/some/place/really/nice'


    >>> expand_template("<install>/really/<type>",\
                         dict(install="/some/place", kind="nice"))
    '/some/place/really/'


    """
    import re
    return re.sub("<(.*?)>", (lambda m: keywords.get(m.group(1),"")), path)

##---------------------------------------------------------------------------##
def convert_jayenne_dependencies(jayenne_deps, keywords):

    """Convert the map of jayenne dependencies into --with-blah=dir
    configure statements. Keywords in the form <fubar> in the path
    specifications are replaced with keywords from the extra names
    parameter arguments.
    """

    strings = []
    for component, path in jayenne_deps.items():

        assert(component in packages)

        expanded_path = expand_template(path, keywords)

        true_path = os.path.normpath(os.path.expanduser(expanded_path))

        is_install_dir(true_path) or sys.exit("Path %s does not appear to \
be a jayenne component installation." % true_path)
        
        strings.append("--with-%s=%s" % (component, true_path))

    return " ".join(strings)

##---------------------------------------------------------------------------##
def _test():
    import doctest
    doctest.testmod()


if __name__=="__main__": 
    _test()

    
