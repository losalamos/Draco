#!/usr/bin/env python

# Python module for interfacing with the module command on ASC
# systems.


import os, Utils

modulepath   = os.environ['MODULEPATH'].split(':')
module_home  = os.environ['MODULESHOME']

def mod_file(s): return (s.find('modulefiles') > -1)
def ver_file(s): return (s.find('version')     > -1)

module_paths    = filter(mod_file , modulepath)
module_versions = filter(ver_file , modulepath)

def module_command(command):

    """Pass the given command to '/usr/bin/modulecmd python' and
    execute the results."""

    exec os.popen('/usr/bin/modulecmd python %s' % command).read()


def avail_modules():

    """Return list of the available modules from the contents of
    os.environ['MODULEPATH']. We only dive one directory level from
    this path in search of files which represent modules.

    """

    modules = []

    for path in module_paths:
        modules.extend(Utils.listdir(path, 1))

    return modules

    
def loaded_modules():
    """Get the list of currently loaded modules."""

    return os.environ['LOADEDMODULES'].split(':')


def add_module(module):
    """Add the module with the provided name"""

    if not module in avail_modules():
        raise "Module %s does not exist." % module
    if not module in loaded_modules():
        module_command("add %s" % module)

    assert(module in loaded_modules())

def remove_module(module):

    if not module in avail_modules():
        raise "Module %s does not exist." % module
    if module in loaded_modules():
        module_command("remove %s" % module)

    assert (module not in loaded_modules())


def list():

    print "Currently Loaded modules:"
    i = 0
    for module in loaded_modules():
        i += 1
        print " %s) %s" % (i, module)

def avail():

    print "Available modules:"
    for module in avail_modules():
        print " %s" % module
    

if __name__=='__main__':

    # Try adding idl. Remove it first to make sure we added it.
    remove_module('tools/idl-6.1')
    list()
    add_module('tools/idl-6.1')
    list()
    avail()

