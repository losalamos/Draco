#!/usr/bin/env python

# Python module for interfacing with the module command on ASC
# systems.


import os, Utils

class ModuleError(Exception): pass


modulepath  = os.environ['MODULEPATH'].split(':')
module_home = os.environ['MODULESHOME']

module_paths = [path for path in os.environ['MODULEPATH'].split(':')
                if (path.find('modulefiles') > -1) and
                os.access(path, os.R_OK)]

avail_modules = reduce( (lambda a,b: a+b),
                        [[module for module in Utils.listdir(path, 1) if
                          '/' in module] for path in
                         module_paths] )

##---------------------------------------------------------------------------##
def module_command(command):

    """Pass the given command to '/usr/bin/modulecmd python' and
    execute the results."""

    exec os.popen('/usr/bin/modulecmd python %s' % command).read()


##---------------------------------------------------------------------------##
def loaded_modules():
    """Get the list of currently loaded modules."""

    try:
        return os.environ['LOADEDMODULES'].split(':')
    except KeyError:
        return []


##---------------------------------------------------------------------------##
def add_module(module):
    """Add the module with the provided name"""

    if not module in avail_modules:
        raise ModuleError ("Module %s does not exist." % module)
    if not module in loaded_modules():
        module_command("add %s" % module)

    assert(module in loaded_modules())

##---------------------------------------------------------------------------##
def remove_module(module):

    if not module in avail_modules:
        raise ModuleError("Module %s does not exist." % module)
    if module in loaded_modules():
        module_command("remove %s" % module)

    assert (module not in loaded_modules())


##---------------------------------------------------------------------------##
def list():

    print "Currently Loaded modules:"
    for i,module in enumerate(loaded_modules()):
        print " %s) %s" % (i+1, module)

##---------------------------------------------------------------------------##
def avail():

    print "Available modules:"
    for module in avail_modules:
        print " %s" % module
    



##---------------------------------------------------------------------------##
def _test():

    # Try adding idl. Remove it first to make sure we added it.
    remove_module('idl/6.1')
    list()
    add_module('idl/6.1')
    list()

if __name__=='__main__':

    _test()
