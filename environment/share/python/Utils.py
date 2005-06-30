#!/usr/bin/env python
#======================================================================
# module: Utils
#
# Contains misc. useful functions
#======================================================================


##---------------------------------------------------------------------------##
def is_attrib(value):
    """is_attrib(value)
    
    Returns true for valid XML attribute specifiers: XXX='YYY'

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
def padlist(l, length):
    """ Pads a list on the right with values of None until it reaches
    the given length

    >>> padlist([1,2,3], 5)
    [1, 2, 3, None, None]
    
    >>> padlist([1,2,3,4,5], 3)
    [1, 2, 3, 4, 5]
    """

    ll = l[:]
    ll.extend([None]*(length-len(l)))
    return ll

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
    regex = re.compile(regexString)

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
def processSetters(options, setters, parameters):

    """
    processSetters(options, setters, parameters)

    "Setters" are command line arguments whose purpose is to provide a
    value for a parameter, i.e. "-x <argument>" means set a
    certain parameter to <argumemt>.

    This function processes a list of options in the form provided by
    getopt. (A list of '(option, value)' tuples). It maps command line
    options to parameters via a dictionary and stores the resulting
    parameter values in another dictionary.

    The argument 'options' contains the list of options from
    getopt.getopt. This consists of a list of tuples where each tuple
    is a pair '(option, value)'.

    Argument 'setters' is the dictionary that maps command line options to
    the parameters that they effect, e.g. {'-x':'x_param'}

    Parameters is the dictionary into which the results of processing
    the command line arguments is written. e.g.
    {'x_param':'value_of_x'} Entries in this argument that are not given
    new values will remain. This is a good mechanism for default
    values.
    

    Examples:

    Map '-x' to x_param and provide it with value: 'x_value':
    >>> params = {}
    >>> processSetters([('-x','x_value')], {'-x':'x_param'}, params)
    >>> print params
    {'x_param': 'x_value'}


    Override the default value for 'x_param': 'x_default'
    >>> params = {'x_param':'x_default'}
    >>> processSetters([('-x','x_value')], {'-x':'x_param'}, params)
    >>> print params
    {'x_param': 'x_value'}
    """

    for option, arg in options:
        for setter, param in setters.items():
            if option == setter: parameters[param] = arg

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
def complete(value, targets):
    """
    Return a list of strings from targets which begin with the
    characters in value
    
    >>> complete('g', ['tall', 'grande', 'venti', 'giant'])
    ['grande', 'giant']
    
    >>> complete('s', ['tall', 'grande', 'venti', 'giant'])
    []
    """

    return [target for target in targets if target.startswith(value)]

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
def _test():
    import doctest, Utils
    doctest.testmod(Utils)

##---------------------------------------------------------------------------##
if __name__=='__main__':
    _test()
