#!/usr/bin/env python
#======================================================================
# module: Utils
#
# Contains misc. useful functions
#======================================================================


##---------------------------------------------------------------------------##
## Function is_attrib:
##---------------------------------------------------------------------------##
def is_attrib(value):
    """Matches only valid XML attribute specifiers: XXX='YYY'

    >>> is_attrib("tag='value'")
    1

    >>> is_attrib('tag="The value is moot!"')
    1

    >>> is_attrib("tag!='value'")
    0

    The tag is limited to characters [A-Za-z] while the value can
    contain anything. Single or double quotes around the tag are okay,
    so long as they match.
    
    """
    import re
    r = re.compile(r'[A-Za-z]+=([\'|\"]).*\1$')

    return bool(r.match(value))


def parse_attrib(value):
    """ Extract the key and value in an XML attribute specifier.

    """
    import re

    if not is_attrib(value):
        raise ValueError, "%s is not a valid XML attribute" % value

    r = re.compile(r'(?P<tag>[A-Za-z]+)=([\'|\"])(?P<value>.*)\2$')

    match = r.match(value)
    return match.group('tag', 'value')


##---------------------------------------------------------------------------##
## Function padList
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
# Function ScanTo
#----------------------------------------------------------------------

def ScanTo(fileHandle, regexString):

    """
    Reads lines in the file pointed to by fileHandle until one matching
    regexString is found. Returns the matching position object and
    leaves the fileHandle in the state after reading that line, if a
    matching line is found. If no matching line is found, the fileHandle
    is left pointing at the line where it was called and the return
    argument is None.
    
    The value of the line can be obtained from the returned position
    object as the .string member.
    """

    import re

    regex = re.compile(regexString)

    # Store the original location for restoration:
    filePosition = fileHandle.tell()

    line = fileHandle.readline()
    while line:
        position = re.search(regexString, line)
        if position:
            break
        else:
            line = fileHandle.readline()

    # If the line is empty, we reached the end of the file without
    # finding the pattern. Before exiting, restore fileHandle to
    # original position in file:

    if line=='': fileHandle.seek(filePosition)

    return position


##---------------------------------------------------------------------------##
def openFileOrString(source):                  

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
class AmbiguousKeyError(Exception):
    """An exception raised by Utils.disambiguate when multiple matching
    values are found."""
    pass

class InvalidKeyError(Exception):
    """An exception raised by Utils.disambiguate when no matching
    values are found."""
    pass

##---------------------------------------------------------------------------##
def disambiguate(value, targets):

    """
    Return a single string from the list argument 'targets' which
    begins with the characters in the string 'value'.

    If no string, or more than one string in targets contains the
    string, raise an "AmbiguousKeyError" exception.

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

def complete(value, targets):
    """
    Return a list of strings from targets which begin with the
    characters in value

   >>> complete('g', ['tall', 'grande', 'venti', 'giant'])
   ['grande', 'giant']
    """

    return [target for target in targets if target.find(value)==0]

##---------------------------------------------------------------------------##
def _test():
    import doctest, Utils
    doctest.testmod(Utils)

##---------------------------------------------------------------------------##
if __name__=='__main__':
    _test()
