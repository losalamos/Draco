import exceptions, os
import Verbosity, Utils, Repo

"""Package CVS

Facilities for interacting with CVS in Python.

Contains classes:

  Repository
  Module
  WorkingCopy

Here's the relationship:

  WorkingCopy ===> Module ---> Repo
              \           \--> Tag (string)
               \     
                \--> Destination


"""

##---------------------------------------------------------------------------##
class ArgumentError(Exception):
    "An exception class for inconsistent combinations of arguments."
    pass

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

class Repository(object):
    """Represents a repository

    A Repository object represents a specific CVS repository.  It
    verifies the existence and readability of the directory where it
    lives.

    >>> r = Repository('/ccs/codes/radtran/cvsroot')
    >>> print r.location
    /ccs/codes/radtran/cvsroot

    >>> print r
    -d /ccs/codes/radtran/cvsroot

    """
    
    def __init__(self, repository):
        self.location = repository
        assert(os.access(self.location, (os.F_OK | os.R_OK)))

    def __str__(self): return "-d %s" % self.location

def lookup_repo(module_name): return Repository(Repo.get_dir(module_name))

##---------------------------------------------------------------------------##
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

class Module(object):
    """Represents a CVS module, versioned with a tag string.

    >>> m = Module('draco/environment', "-r dummy_tag")
    >>> m.name
    'draco'
    >>> m.module
    'draco/environment'
    >>> print m
    -r dummy_tag draco/environment
    
    """
    
    def __init__(self, module_name, tag):
        self.tag     = tag
        self.module  = module_name
        self.name    = module_name.split("/",1)[0]

        # Lookup the correct repository:
        self.repo = lookup_repo(self.name)

    def __str__(self): 
        return "%s %s" % (self.tag, self.module)


##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

class WorkingCopy(Module):

    """A WorkingCopy is a module, assigned to a particular directory
    name for checking out.

    >>> w = WorkingCopy('draco/environment', "-r dummy_tag", 'environment')
    >>> w.output_dir()
    'environment'
    >>> w.checked_out()
    False

    For now, all you can do with it is check it out.
    """

    def __init__(self, module_name, tag, destination=None):

        super(WorkingCopy, self).__init__(module_name, tag)

        self.destination = destination
        self.path        = None

    def __str__(self):
        if self.destination:
            part = "-d %s " % self.destination
        else:
            part = " "

        return part + super(WorkingCopy,self).__str__()
        
    def output_dir(self):  return self.destination or self.name

    def checked_out(self): return bool(self.path)

    def checkout(self, location, export=False, verbose=Verbosity.ignore()):
        # Switch to the indicated directory.
        try:
            os.chdir(location)
        except OSError:
            sys.exit("Could not chdir to directory %s" % location)
            
        cvs_command = export and "export" or "checkout"

        command = "cvs -Q %s %s %s" % (self.repo, cvs_command, self)

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


def make_working_copy(module_name, checkout_name, tag_kind, tag_name):
    tag = make_tag(tag_kind, tag_name)
    return WorkingCopy(module_name, tag, checkout_name)


##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
def _test():
    import doctest, CVS
    return doctest.testmod(CVS)


if __name__=="__main__":
    _test()
