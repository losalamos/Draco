import Jayenne, Utils

""" These are functions common to applications that work with CVS
repositories. Currently: checkout and pitch.

"""


##---------------------------------------------------------------------------##
def parse_input(options, args):
    """Build a CVS.Module object from the options and arguments

    Supported Module/tag combos:

    1.  module              | module is args[0], tag is not specified -> HEAD.
    2.  module -r tag       | tag is options.tagname, module is args[0]
    2b. module/part -r tag  | same as above, but module is args[0].split('/')[0]
                            | CVS.Module takes care of this case.
    3.  -r module-version   | Tag is options.tagname, module is
                            | options.tagname.split('-')[0]
    4.  module version      | Both in args. Tag is 'module-version'
    4b. module/part version | Tag is built from 'module' not 'module/part'
    
    """

    if len(args) > 2:
        raise Utils.ArgumentsException("Too many arguments.")

    
    # If there is nothing in agrs, we must be in case 3.
    if not args:
        parts = options.tagname.split("-",1)
        if len(parts) < 2:
            raise Utils.ArgumentsException( \
                "Missing package name. Could not deduce from tag.")
        else:
            module = parts[0]
            tag    = options.tagname
            
    else:
        module = args[0]

        if options.tagname:
            tag = options.tagname
        elif len(args)==2:
            tag = module.split('/',1)[0] + '-' + args[1]
        else:
            tag=None

    
    return (module, tag)


##---------------------------------------------------------------------------##
def parse_tag(tag, splitter=":"):
    """Extract the kind and label of a tag from the text
    description"""

    if not tag: return ('head', None)

    (prefix, name) = Utils.padlist(tag.split(splitter, 1), 2)

    if not name: return ('symbolic', prefix)

    return (prefix, name)
    

##---------------------------------------------------------------------------##
def checkout(module, dir_name, path_name, export, tag_kind, tag_name, verbosity):
    "Get a WorkingCopy object and check it out"

    work_copy = Jayenne.make_working_copy(module, dir_name, tag_kind, tag_name)
    work_copy.checkout(path_name, export, verbosity)

    return work_copy


