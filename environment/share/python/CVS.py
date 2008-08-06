import exceptions, os
import Verbosity, Utils, Repo

"""Package CVS

Facilities for interacting with CVS in Python.

Contains classes:

  Repository
  Package
  WorkingCopy

Here's the relationship:

  WorkingCopy ===> Package ---> Repo.Repositiry
              \           \---> Tag (string)
               \     
                \--> Destination


"""

##---------------------------------------------------------------------------##
class ArgumentError(Exception):
    "An exception class for inconsistent combinations of arguments."
    pass

##---------------------------------------------------------------------------##

def make_tag(kind, name=None):
    """Convert a tag kind and name into a tag.  Performs completion on
    the kind.

    This could be a static member of the Tag class.
    """
    kinds = ['head', 'date', 'symbolic']

    try:
        kind = Utils.disambiguate(kind, kinds)
    except Utils.KeyError, e:
        raise ArgumentError("Bad tag prefix %s" % e.args[0])

    if not name and kind != 'head':
        raise ArgumentError("Date and revision tags need a name")

    if name and kind == 'head':
        raise ArgumentError("Head tags to net have a name")

    if   kind=='head':     return ""
    elif kind=='date':     return "-D %s" % name
    elif kind=='symbolic': return "-r %s" % name
    else:
        raise ArgumentError("Unrecognized tag kind: %s" % kind) 
    
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

class Package(object):
    """Represents a code package, defined as a subdirectory of the CVS
    root directory continaing source code.

    Because a package is contained in a single directory in the CVS
    directory tree, we can copy it to populate another repository.

    This is less general than a CVS module, which is anything that can
    be checked out, and comes in several kinds.

    >>> m = Package('draco/environment', Utils.Bunch(location='stub'))
    >>> m.name
    'draco/environment'
    >>> m.path
    'stub/draco/environment'
    >>> print m
    draco/environment
    
    """
    
    def __init__(self, package_name, repository):
        self.name       = package_name
        self.repository = repository
        self.path       = os.path.join(self.repository.location, self.name)
        
    def __str__(self): 
        return self.name
    

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

class Module(object):
    """Represents a CVS module. This is basically anything that can be
    checked out of CVS and is usually defined in the CVSROOT/module
    file.

    Have to: 

    Figure out how to get the reposotiry from the module name. Put
    these module names in the master table also?

    Verify the validity of the module name.
    """
    pass


##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

class WorkingCopy(object):

    """A WorkingCopy is a package, with an optional versioning tag,
    and assigned to a location for checking out, either at
    construction or when checked out.

    >>> w = WorkingCopy('draco/environment', "-r dummy_tag", 'environment')
    >>> w.output_dir()
    'environment'
    >>> w.checked_out()
    False

    For now, all you can do with it is check it out.

    """

    def __init__(self, package, tag, destination=None):

        self.package     = package
        self.tag         = tag
        self.destination = destination
        self.path        = None

    def __str__(self):
        if self.destination:
            part = "-d %s " % self.destination
        else:
            part = " "

        return "%s %s %s" % (part, self.tag, self.package)
        
    def output_dir(self):  return self.destination or self.package.name

    def checked_out(self): return bool(self.path)

    def checkout(self, location, export=False, verbose=Verbosity.ignore()):

        try:
            os.chdir(location)
        except OSError:
            sys.exit("Could not chdir to directory %s" % location)
            
        cvs_command = export and "export" or "checkout"

        command = "cvs -Q %s %s %s" % (
            self.package.repository, 
            cvs_command, 
            self.__str__())


        verbose("Executing CVS command: %s" % command, 1)
        verbose("in directory %s" % location, 2)

        command_out = os.popen(command)
        output = command_out.read()
        error_code = command_out.close()
    
        if error_code:
            raise exceptions.RuntimeError(
                "CVS command failed with error: %s" % error_code)

        self.path = os.path.join(location, self.output_dir());
        assert(os.access(self.path, os.F_OK | os.R_OK))
        
        return self.path


##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
def _test():
    import doctest, CVS
    return doctest.testmod(CVS)


if __name__=="__main__":
    _test()
