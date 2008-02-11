"""A module for setting up the configuration and building of Jayenne
codes.
""" 

class PackageError(Exception):
    """Reports errors relating to the Jayenne packages"""
    pass

class LookupError(PackageError): pass


# Ugh. This information is repeated in Repo.
packages     = ['draco', 'clubimc', 'milagro', 'wedgehog']

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
    LookupError: Found no packages
    """

    results = [package for (package,comp_list) in components.items() if
               component in comp_list]

    if len(results) > 1: 
        raise LookupError("Found more than one package")
    if len(results) < 1:
        raise LookupError("Found no packages")
    
    return results[0]


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
def make_dependency_string(component, path):
    """Return a configure string which specified a jayenne component
    dependency. """

    assert(component in packages)

    return "--with-%s=%s" % (component, path)


def _test():
    import doctest
    doctest.testmod()


if __name__=="__main__": 
    _test()

    
