import os
import fnmatch

"""This module contains small functions from David Beazley's Talk
"Generator Tricks for Systems Programmers" at PyCon'2008 and some
functions based on the same ideas.

See http://www.dabeaz.com/generators/ for the complete talk and source
code. 
"""


##---------------------------------------------------------------------------##
def gen_find(filepat, top, recurse=True):
    """Generate (yield) files in directory tree under top which satisfy
    filepat.
    """

    for path, dirlist, filelist in os.walk(top):
        if not recurse: dirlist[:]=[]
        for name in fnmatch.filter(filelist,filepat):
            yield os.path.join(path,name)


##---------------------------------------------------------------------------##
import gzip, bz2

def gen_open(filenames):
    """Yeild a filehandle for the filenames in the argument. Use gzip
    or bz2 modules as needed for opening.
    """

    for name in filenames:
        if name.endswith(".gz"):
            yield gzip.open(name)
        elif name.endswith(".bz2"):
            yield bz2.BZ2File(name)
        else:
            yield open(name)

##---------------------------------------------------------------------------##
import re

def gen_grep(pat,lines):
    """Yeild all of the lines in argument lines which match regular
    expression in pat.
    """
    patc = re.compile(pat)
    for line in lines:
        if patc.search(line): yield line


##---------------------------------------------------------------------------##
def gen_subdirs(top, recurse=True):
    """Generate the sub-directories under top"""

    import os.path
    
    for path, dirlist, filelist in os.walk(top):
        for name in dirlist:
            yield os.path.join(path,name)
        if not recurse: dirlist[:]=[]

