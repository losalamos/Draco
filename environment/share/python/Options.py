#!/usr/bin/env python
#======================================================================
# module: Options
#
# Contains misc. useful functions for processing options
#======================================================================

import Utils

##---------------------------------------------------------------------------##
class Exclusion(object):

    """A helper class to keep a list of options which may not appear
    together and check/manipulate a bunch object into compliance.
    
    >>> e = Exclusion(['opt','debug'], 'debug')
    >>> e.options
    ['opt', 'debug']
    >>> e.default_option
    'debug'

    An exclusion can involve a value option:
    >>> e = Exclusion(['verb:', 'quiet'], 'quiet')
    >>> p = Utils.Bunch()
    >>> e.check_params(p)
    >>> hasattr(p, 'quiet')
    1
    
    The default can also hold a value:
    >>> e = Exclusion(['cache:', 'no_cache'], 'cache:10')
    >>> e.default_option
    'cache:'
    >>> e.default_value
    '10'
    >>> p = Utils.Bunch()
    >>> e.check_params(p)
    >>> int(p.cache)
    10

    >>> p.no_cache = None
    >>> e.check_params(p)
    Traceback (most recent call last):
    ...
    ValueError: Mutually exclusive options "cache:" and "no_cache".
    """

    def __init__(self, options, default):

        self.options = options

        self.default_option, self.default_value = \
                             Utils.padlist(default.split(':'),2)

        if self.default_value: self.default_option += ":"

        if self.default_option not in self.options:
            raise ValueError, "Default not in options"


    # Strip off the colon of attrib if necessary
    def _attribute_name(self, attrib): return attrib.rstrip(':')

    def _set_default(self, params):
        setattr(params, self._attribute_name(self.default_option), self.default_value)
        
    def check_params(self, params):

        """Ensure that no more than one of the mutaully exclusive
        options appear. If none appear, set the default.

        """

        # Count the number of self.options which appear as attributes
        # of params.
        opt_matches = []
        for option in self.options:
            if hasattr(params, self._attribute_name(option)): opt_matches.append(option)
            
            # If more than one, raise an error
            if len(opt_matches) > 1:
                raise ValueError, "Mutually exclusive options \"%s\" " \
                      "and \"%s\"." % (opt_matches[0], opt_matches[1]) 

        # If none, set the default
        if len(opt_matches) == 0: self._set_default(params)

            
##---------------------------------------------------------------------------##
def process_options(options):

    """Convert a list of options into a bunch object.

    Options can be simple: 'name', or a labeled value: 'name:value'.

    A simple option is stored in the bunch object with a correspinding
    value of None.  The value of a simple option stored in a parameter
    object does not matter. It should only be checked when it makes
    sense to do so.

    >>> options = process_options(['cache:10', 'opt'])
    >>> options.cache
    '10'
    >>> hasattr(options,'opt')
    1
    
    """

    d =dict( [Utils.padlist(option.split(':'),2) for option in
              options] )
    

    return Utils.Bunch(d)



##---------------------------------------------------------------------------##
def extract_value(arguments, label, default = None, type = lambda x:x):

    """ Extracts the value from the matching label:value pair in a
    list 'args'

    Multiple occurances cause the last matching value to be
    returned.

    No occurances cause None to be returned, unless optional argument
    'default' is provided.

    Optional argument 'type' can be any callable object and it applied
    to the value before returning it, if it is found. Otherwise it is
    appled to the default argument, if provided. 

    >>> args = Utils.Bunch(n='10', label='value')
    >>> args.n
    '10'
    >>> args.label
    'value'

    >>> extract_value(args, 'label')
    'value'

    >>> extract_value(args, label='n', type = int)
    10

    >>> print extract_value(args, 'wrong', type = lambda x: x**2)
    None

    >>> print extract_value(args, 'wrong', default = 10)
    10

    """

    try:
        return type(getattr(arguments, label))
    except AttributeError:
        if default:
            return type(default)
        else:
            return None
        
    
##---------------------------------------------------------------------------##
def _test():
    import doctest, Options
    doctest.testmod(Options)

##---------------------------------------------------------------------------##
if __name__=='__main__':
    _test()
