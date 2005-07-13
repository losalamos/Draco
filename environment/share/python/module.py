#!/usr/bin/env python

# Python module for interfacing with the module command on ASC
# systems.


import os, Utils

modulepath = os.environ['MODULEPATH'].split(':')
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
    """Get the list of modules currently loaded."""

    return os.environ['LOADEDMODULES'].split(':')
    

if __name__=='__main__':

    # Try adding tool. Remove it first to make sure we added it.
    module_command('remove tools/idl-6.1')
    module_command('list')
    module_command('add tools/idl-6.1')

    print loaded_modules()

    print avail_modules()


