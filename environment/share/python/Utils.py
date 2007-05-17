#!/usr/bin/env python
#======================================================================
# module: Utils
#
# Contains misc. useful functions
#======================================================================

##---------------------------------------------------------------------------##
## Basic data types:
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
class Bunch(object):
    """An object created with attributes from a dictionary and/or a
    keyword list

    >>> d = {'key': 'value'}
    >>> b = Bunch(d, second_key='another_value')
    >>> b.key
    'value'
    >>> b.second_key
    'another_value'

    """

    def __init__(self, dict = {}, **kwds):
        self.__dict__.update(dict)
        self.__dict__.update(kwds)



##---------------------------------------------------------------------------##
class Parameters(object):
    """A specialized kind of Bunch which keeps a list of expected
    parameters and can print their values.

    >>> params = ['key', 'second_key']
    >>> d = {'key': 'value'}
    >>> b = Parameters(params, d, second_key='another_value')
    >>> print b
             key: value
      second_key: another_value

    Ideally, this class would be derived from Bunch, but I was not
    able to figure out how to pass the key-word arguments to the base
    class constructor.

    """

    def __init__(self, params, dict = {}, **kwds):
        self.params = params
        self.__dict__.update(dict)
        self.__dict__.update(kwds)

    def __str__(self):
        s = ''
        display_width = max(map(len, self.params)) + 2

        for param in self.params:
            s += "%*s: %s" % (display_width, param, getattr(self, param))
            s += '\n'

        return s[:-1]


##---------------------------------------------------------------------------##
class Data(object):
    """ Accumulate and report statistics on data.

    >>> d = Data()
    >>> d.add(10)
    >>> d.add(20)
    >>> d.add(30)
    >>> d.trials
    3
    >>> d.average()
    20.0
    >>> d.deviation()
    10.0
    """
    
    def __init__(self):
        self.values = []
        self.value  = 0.0
        self.trials = 0

    def add(self, value):
        self.values.append(value)
        self.value  += value
        self.trials += 1

    def average(self):
        return self.value / self.trials

    def deviation(self):

        if self.trials < 2: return 0.0

        average = self.average()

        import math

        deviations = map(lambda x: (x-average)**2, self.values)
        summation  = reduce(lambda x,y: x+y, deviations)

        return math.sqrt(summation / ( self.trials - 1))



##---------------------------------------------------------------------------##
## List manipulations
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
def padlist(l, length, value=None):

    """ Pads a list on the right with 'value' (default None) until it
    reaches the given length

    >>> padlist([1,2,3], 5)
    [1, 2, 3, None, None]
    
    >>> padlist([1,2,3], 5, 0)
    [1, 2, 3, 0, 0]
    
    >>> padlist([1,2,3,4,5], 3)
    [1, 2, 3, 4, 5]
    """

    ll = l[:]
    ll.extend([value]*(length-len(l)))
    return ll

##---------------------------------------------------------------------------##
def unique_append(a_list, an_item):

    """Append an_item to a_list only if it does not yet appear

    >>> unique_append([1,2,3],4)
    [1, 2, 3, 4]

    >>> unique_append([1,2,3,4],4)
    [1, 2, 3, 4]
    """
    if an_item not in a_list: a_list.append(an_item)

    return a_list

##---------------------------------------------------------------------------##
def unique_extend(a_list, b_list):

    """Extend a_list with items in b_list if they do not already
    appear.

    >>> unique_extend([1,2,3], [3,4,5])
    [1, 2, 3, 4, 5]

    """

    for an_item in b_list:
        if an_item not in a_list: a_list.append(an_item)

    return a_list




##---------------------------------------------------------------------------##
## XML stuff:
##---------------------------------------------------------------------------##



##---------------------------------------------------------------------------##
def is_attrib(value):
    """is_attrib(value)
    
    Returns true for valid XML attribute specifiers

    >>> is_attrib("tag='value'")
    1

    >>> is_attrib('tag="value"')
    1

    >>> is_attrib("tag,'value'")
    0

    The tag is limited to characters [A-Za-z] while the value can
    contain anything. Single or double quotes around the tag are okay,
    so long as they match.
    
    """
    import re
    r = re.compile(r'[A-Za-z]+=([\'|\"]).*\1$')

    return bool(r.match(value))


##---------------------------------------------------------------------------##
def parse_attrib(value):
    """parse_attrib(value)

    Returns as a pair the key and value in an XML attribute specifier.

    >>> parse_attrib("tag='value'")
    ('tag', 'value')

    >>> parse_attrib('tag="value"')
    ('tag', 'value')

    >>> parse_attrib("tag,value")
    Traceback (most recent call last):
    ...
    ValueError: 'tag,value' is not a valid XML attribute.

    """
    import re

    r = re.compile(r'(?P<tag>[A-Za-z]+)=([\'|\"])(?P<value>.*)\2$')

    match = r.match(value)
    if not match:
        raise ValueError, "'%s' is not a valid XML attribute." % value
    else:
        return match.group('tag', 'value')


##---------------------------------------------------------------------------##
## File manipulation:
##---------------------------------------------------------------------------##

#----------------------------------------------------------------------
def ScanTo(fileHandle, regexString):

    """
    ScanTo(fileHandle, regexString)

    Reads lines in the file pointed to by fileHandle until one
    matching regexString is found. If a match is found, returns the
    match object and leaves the fileHandle at pointing to the next
    line. If no matching line is found, the fileHandle is left
    pointing at the line where it was called and the return argument
    is None.
    
    In this example, we scan this file to find the beginning of this
    function, then scan for the import command and retrieve the module
    name. We then fail at another scan and verify that the filehandle
    is in the right place.

    This test must be run from the directory containing Utils.py

    >>> file = open('Utils.py')
    >>> match = ScanTo(file, '^[ ]*def[ ]+ScanTo')
    >>> match.string.strip()
    'def ScanTo(fileHandle, regexString):'

    >>> match = ScanTo(file, '^[ ]*import[ ]+(.*)')
    >>> match.group(1)
    're'

    >>> match = ScanTo(file, "^Can\'t match this.$")
    >>> print match
    None
    >>> file.readline().strip()
    'regex = re.compile(regexString)'

    """

    import re
    if isinstance(regexString, str):
        regex = re.compile(regexString)
    else:
        regex = regexString

    original_position = fileHandle.tell()
    
    match = None
    while not match:
        line = fileHandle.readline()
        match = regex.search(line)
        if not line: break

    if not match: fileHandle.seek(original_position)

    return match


##---------------------------------------------------------------------------##
def openFileOrString(source):
    """openFileOrString(source)

    A unified open function which will open a specified file, if it
    exists, or, failing that, treat the argument as a string and open
    it for reading.

    >>> file = openFileOrString('Utils.py')
    >>> file.readline().strip()
    '#!/usr/bin/env python'

    >>> file = openFileOrString('This is a test.')
    >>> file.readline().strip()
    'This is a test.'

    """

    # try to open with native open function (if source is pathname)
    try:                                  
        return open(source)                
    except (IOError, OSError):            
        pass                              
    
    # treat source as string
    import StringIO                       
    return StringIO.StringIO(str(source))  


##---------------------------------------------------------------------------##
## Option processing
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
def process_setter_options(options, option_param_map, defaults):

    """
    process_setter_options

    Processes command-line option pairs: ('-option', 'value'), by storing
    the values as data members of an object.

    Argument options is a list of the command-line options pairs taken
    from the command line like, e.g.:

    >>> options = [('-x', 'x_value'), ('-a', 'a_value')]

    Argument setter_map is a dictionary which maps command-line
    options to the parameters that they set. e.g:

    >>> option_param_map = {'-x': 'x_param', '-y': 'y_param', '-z': 'z_param'}

    Argument defaults is a dictionary which maps parameters to default
    values. Parameters without default values default to None

    >>> defaults = {'x_param': 'x_default', 'y_param': 'y_default'}

    The return value is a 'bunch' containing the variables as members
    with the assigned values and a list of the unrecognized options. e.g:

    >>> (params, extras) = process_setter_options(options, option_param_map, defaults)

    >>> print params.x_param
    x_value

    >>> print params.y_param
    y_default

    >>> print params.z_param
    None

    >>> extras
    [('-a', 'a_value')]

    """

    # Extract the list of expected parameters:
    parameters = option_param_map.values()

    # Add a default of None for all parameters not already in defaults
    for parameter in parameters:
        defaults.setdefault(parameter, None)

    # A list of options we did not recognize:
    unknowns = []
        
    # Loop over the (option, value) pairs in options. Lookup the
    # corresponding parameters from option_param_map. Modify parameter
    # values in place in the default dictionary.

    for option, value in options:
        parameter = option_param_map.get(option)
        if parameter:
            defaults[parameter] = value
        else:
            unknowns.append((option, value))
        
            

    # Convert the default dictionary into a Parameters object. Return it and the
    # unknown arguments.

    return Parameters(parameters, defaults), unknowns





##---------------------------------------------------------------------------##
## String completion
##---------------------------------------------------------------------------##


##---------------------------------------------------------------------------##
class KeyError(Exception):
    "An exception class for all disambuguation errors."
    pass

class AmbiguousKeyError(KeyError):
    """An exception raised by Utils.disambiguate when multiple matching
    values are found."""
    pass

class InvalidKeyError(KeyError):
    """An exception raised by Utils.disambiguate when no matching
    values are found."""
    pass

##---------------------------------------------------------------------------##
def complete(value, targets):
    """
    Return a list of strings from targets which begin with the
    characters in value
    
    >>> complete('g', ['tall', 'grande', 'venti', 'giant'])
    ['grande', 'giant']
    
    >>> complete('gr', ['tall', 'grande', 'venti', 'giant'])
    ['grande']
    
    >>> complete('s', ['tall', 'grande', 'venti', 'giant'])
    []
    """

    return [target for target in targets if target.startswith(value)]

##---------------------------------------------------------------------------##
def disambiguate(value, targets):

    """
    Return a single string from the list argument 'targets' which
    begins with the characters in the string 'value'.

    If more than one string in targets contains the string, raise an
    "AmbiguousKeyError" exception. If no strings match, raise
    "InvalidKeyError".

    >>> disambiguate('gr', ['tall', 'grande', 'venti', 'giant'])
    'grande'

    Ambiguous values generate an exception:
    >>> disambiguate('g', ['tall', 'grande', 'venti', 'giant'])
    Traceback (most recent call last):
    ...
    AmbiguousKeyError: ('g', ['tall', 'grande', 'venti', 'giant'])

    Failure to match any any generates an exception:
    >>> disambiguate('medium', ['tall', 'grande', 'venti', 'giant']) 
    Traceback (most recent call last):
    ...
    InvalidKeyError: ('medium', ['tall', 'grande', 'venti', 'giant'])
    """

    matches = complete(value,targets)

    if len(matches)==1:
        return matches[0]
    elif len(matches) > 1:
        raise AmbiguousKeyError(value, targets)
    else:
        raise InvalidKeyError(value, targets)





##---------------------------------------------------------------------------##
## Filesystem:
##---------------------------------------------------------------------------##


def get_timestamp_as_string(filename):
    """
    Return the timestamp of the given filename formatted as a string

    >>> get_timestamp_as_string("/home/mwbuksas")
    'Tue Sep 12 15:31:51 2006'
    """

    import os, stat, time
    return time.asctime(time.localtime(os.stat(filename)[stat.ST_MTIME]))


##---------------------------------------------------------------------------##
def listdir(dir, recurse = 0):

    """Get relative names of files in 'dir' up to the specified
    recursion depth. Use a negative recurse argment for arbitrary
    depth.  """

    import os

    names = [os.path.join(dir, name) for name in os.listdir(dir)]

    files = map(os.path.basename, filter(os.path.isfile, names))
    dirs  = filter(os.path.isdir,  names)

    if not recurse == 0:
        for subdir in dirs:

            subdir_name = os.path.basename(subdir)

            sub_files = listdir(subdir, recurse-1)
            sub_files = [os.path.join(subdir_name, name) for name in
                         sub_files]
            files.extend(sub_files)

    return files




##---------------------------------------------------------------------------##
def listFiles(root, patterns='*', recurse=True, relative_paths=False,
              return_folders=False):
    """ From the Python Cookbook version 1 section 4.18 by Robin
    Parmar and Alex Martelli

    Modified by Mike Buksas to add relative_paths option.
    """

    import os.path, fnmatch
 
    # Expand patterns from semicolon-separated string to list
    pattern_list = patterns.split(';')

    # We turn the arguments into a bunch for passing into os.path.walk:
    arg = Bunch(root=root,
                recurse=recurse,
                pattern_list=pattern_list,
                return_folders=return_folders,
                results=[])
 
    # This is the function we pass to os.path.walk. 
    def visit(arg, dirname, files):
        
        # Append to arg.results all relevant files (and perhaps folders)
        for name in files:

            fullname = os.path.normpath(os.path.join(dirname, name))

            if relative_paths:
                display_name = name
            else:
                display_name = fullname

            if arg.return_folders or os.path.isfile(fullname):
                for pattern in arg.pattern_list:
                    if fnmatch.fnmatch(name, pattern):
                        arg.results.append(display_name)
                        break

        # If not recursing, clobber the file list to stop recursion.
        if not arg.recurse: files[:]=[]
 
    os.path.walk(root, visit, arg)
 
    return arg.results


##---------------------------------------------------------------------------##
def _test():
    import doctest, Utils
    doctest.testmod(Utils)

##---------------------------------------------------------------------------##
if __name__=='__main__':
    _test()
