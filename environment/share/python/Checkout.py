import CVS, Utils


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

    work_copy = CVS.make_working_copy(module, dir_name, tag_kind, tag_name)
    verbosity("Identified module %s" % module, 2)
    
    work_copy.checkout(path_name, export, verbosity)

    return work_copy


