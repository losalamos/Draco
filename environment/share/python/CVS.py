import exceptions, os
import Verbosity, Utils

"""Package CVS

Facilities for interacting with CVS in Python.

A Repository objects represents a specific CVS repository. It
implements commands to check out modules from it.
"""

##---------------------------------------------------------------------------##
class ArgumentError(Exception):
    "An exception class for inconsistent combinations of arguments."
    pass


##---------------------------------------------------------------------------##
class Tag:
    """Represents the kinds of tags that we can apply to a package
    when checking it out.

    Converts itself to a string for composition of CVS commands

    >>> print Tag('r', 'tag_name')
    -r tag_name


    >>> print Tag('x', 'tag_name')
    Traceback (most recent call last):
    ...
    ArgumentError: Unrecognized tag kind: x


    >>> t = Tag('D', 'yesterday')
    >>> t.kind
    'D'
    >>> t.name
    'yesterday'

    >>> str(Tag())
    '-r HEAD'

    >>> t = Tag.build('some_branch')
    >>> t.kind
    'r'

    >>> t = Tag.build('date:sometime')
    >>> t.kind
    'D'
    >>> t.name
    'sometime'

    """

    kinds = ['h', 'D', 'r']

    # Add a head kind?
    descriptors = {'date':'D','symbolic':'r'}

    def __init__(self, kind="-r", name=""):

        self.kind = kind
        self.name = name

        # If no name provided, no tag is intended.
        if not self.name:
            self.kind = 'h'

        if (self.name and not self.kind):
                raise ArgumentError("Must provide kind of tag.");

        if self.kind not in self.kinds:
            raise ArgumentError(
                "Unrecognized tag kind: %s" % self.kind)

    def __str__(self): return self.tag_string()

    def tag_string(self):

        if not self.kind:      return "-r HEAD"
        elif self.kind == 'D': return "-D %s" % self.name
        elif self.kind == 'r': return "-r %s" % self.name
        else:                  return "-r HEAD"
        

    def build(desc, splitter=":"):
        "Look for keywords in the tag string which indicate the type."

        parts = desc.split(splitter, 1)

        # If splitter not encountered, assume symbolic tag
        if len(parts) == 1:
            return Tag('r', desc)
        else:
            kind = Utils.complete(parts[0], Tag.descriptors.keys())
            # If the first word fits exactly one key, use it:
            if len(kind)==1:
                return Tag(Tag.descriptors[kind[0]], parts[1])
            else:
                return desc
        
    build = staticmethod(build)            



##---------------------------------------------------------------------------##
class Module:
    """Represents a module in a CVS repository.

    The directory varaible is used to select a destination name for
    the result when checked out.

    >>> m = Module('draco/environment', 'environment')
    >>> m.name
    'draco'
    >>> m.module
    'draco/environment'
    >>> m.destination
    'environment'
    """
    
    def __init__(self, module, destination="", tag="" ):

        self.name        = module.split("/",1)[0]
        self.module      = module
        self.destination = destination
        self.tag         = Tag.build(tag)

    def dir_string(self):

        if self.destination:
            return "-d %s" % self.destination
        else:
            return ""

    # TODO Make this spew more stuff:
    def __str__(self):
        return self.module

        


##---------------------------------------------------------------------------##
class Repository:
    """Represents a repository

    >>> r = Repository('/codes/radtran/cvsroot')
    >>> print r.location
    /codes/radtran/cvsroot

    """
    

    def __init__(self, repository):

        self.location = repository
        assert(os.access(self.location, (os.F_OK | os.R_OK)))




##---------------------------------------------------------------------------##
def checkout(repository,
             module,
             location,
             export,
             verbose = Verbosity.ignore()
             ):
    """Function checkout implemenets the checkout and export cvs
    commands given a description of the desired module, the repository
    it comes from, and where the result should go.
    """

    cvs_command = export and "export" or "checkout"

    command = "cvs -Q -d %s %s %s %s %s" % \
              (repository.location,
               cvs_command,
               module.dir_string(),
               module.tag.tag_string(),
               module.module
               )

    verbose("Executing CVS command: %s" % command, 1)
    verbose("in directory %s" % location, 2)

    
    # Switch to the indicated directory.
    try:
        os.chdir(location)
    except OSError:
        sys.exit("Could not find directory %s" % location)
            
    command_out = os.popen(command)
    output = command_out.read()
    error_code = command_out.close()
    
    if error_code:
        raise exceptions.RuntimeError(
            "CVS command failed with error: %s" % error_code)

    # Figure out what directory we just created. If we used a module
    # destination, return that appended to the location. Otherwise,
    # append the module path to the location and return that

    new_path = module.destination or module.module
    path     = os.path.join(location, new_path);

    assert(os.access(path, os.F_OK | os.R_OK))

    return path


def _test():
    import doctest, CVS
    return doctest.testmod(CVS)


if __name__=="__main__":
    _test()
