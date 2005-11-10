#!/usr/bin/env python
""" Module Select

Contains functions for selecting an object from a collection based on
values of it's attributes
"""

__author__ = "Michael Buksas"
__date__   = "9 November 2005"

##---------------------------------------------------------------------------##
def evaluator(obj, str, object_name="obj", local_vars={}, global_vars={}): 
    """Evaluate a string tailored to extract information about a
    specific object.

    Argument obj is the object we are interested in.

    Argument str is a string to be evaluated.

    Argument object_name is the name of the object as it appears in
    argument str. It defaults to 'obj'

    Arguments locals and globals are other variables to add to the
    local and global environments for the evaluation.

    >>> from Utils import Bunch
    >>> obj = Bunch(name='blah', data=42.6)
    >>> evaluator(obj, "obj.data > 10.0")
    1

    >>> widget = obj
    >>> evaluator(obj, "widget.data > 10.0", "widget")
    1

    >>> evaluator(obj, "obj.name")
    'blah'

    >>> thresh_value = 10.0
    >>> evaluator(obj, "obj.data > thresh", local_vars = {'thresh':thresh_value})
    1


    >>> import math
    >>> evaluator(obj, "math.sin(obj.data)", global_vars = globals())
    -0.98228657290363453

    Or, extract the specific module from the globals
    >>> evaluator(obj, "math.sin(obj.data)", global_vars = {'math':globals()['math']})
    -0.98228657290363453

    >>> evaluator(obj, "math.sin(obj.data)")
    Traceback (most recent call last):
    ...
    NameError: name 'math' is not defined

    """

    local_vars[object_name]=obj

    return eval(str, global_vars, local_vars)


##---------------------------------------------------------------------------##
def has_attrs(obj, attribs):
    """Return true/false if the given object's attributes match the
    values in the attribs dictionary.

    >>> from Utils import Bunch
    >>> o = Bunch(a=1, b=2)

    >>> attribs = {'a':1}
    >>> has_attrs(o, attribs)
    1

    >>> wrong_attrib = {'something': 10}
    >>> has_attrs(o, wrong_attrib)
    0

    >>> wrong_value = {'a': 10}
    >>> has_attrs(o, wrong_value)
    0

    """

    for (key, value) in attribs.items():
        if not hasattr(obj, key) or not (getattr(obj, key) == value):
            return False

    return True

##---------------------------------------------------------------------------##
def satisfies_all(obj, descrs, name="obj", local_vars={}, global_vars={}):
    """Determine if an object satisfies a list of descriptions.

    Argument obj is the object to test

    Argument descrs is a list of strings containing evaluatable python
    expressions. These will be applied to the object in argument obj

    Argument name is the name of the object as is appears in the
    strings in argument descr.

    Arguments locals and globals are variable dictionaries which can
    contain other information necessary to evaluate the expressions in
    argument descr.

    >>> epsilon = 0.1
    >>> name = "the one"

    >>> descrsA = ['obj.data > epsilon']
    >>> descrsB = ['obj.name == name', 'obj.data > 0']

    >>> from Utils import Bunch
    >>> obj1 = Bunch(data=10.0, name="some object")
    >>> obj2 = Bunch(data=0.01, name="the one")

    >>> variables = {'epsilon':epsilon, 'name': name}

    >>> satisfies_all(obj1, descrsA, local_vars = locals())
    1

    >>> satisfies_all(obj1, descrsB, 'obj', variables)
    0
    
    >>> satisfies_all(obj2, descrsB, local_vars = variables)
    1

    """

    for descr in descrs:
        if not evaluator(obj, descr, name, local_vars, global_vars):
            return False

    return True

##---------------------------------------------------------------------------##
def exact_filter(objects, attribs):
    """Return a new list of objects extracted from argument
    'objects' which satisify a specific set of attribute requirements.

    Argument objects contains the list of objects

    Argument arrtibs is a dictionary of attribute, value pairs

    Objects are selected for the new list if all of their attributes
    which are keys in the attribs argument match the corresponding
    values.

    >>> from Utils import Bunch
    >>> o1 = Bunch(a=1, b=2)
    >>> o2 = Bunch(a=1, c=3)
    >>> o = [o1,o2]

    >>> attribs_both    = {'a':1}
    >>> attribs_1       = {'b':2}
    >>> attribs_neither = {'b':2, 'c':3}

    >>> len(exact_filter(o, attribs_both))
    2

    >>> len(exact_filter(o, attribs_1))
    1

    >>> len(exact_filter(o, attribs_neither))
    0
    
    """

    return [object for object in objects if has_attrs(object, attribs)]


##---------------------------------------------------------------------------##
def _test():
    import doctest, Select
    doctest.testmod(Select)



if __name__=='__main__':
    _test()
    
