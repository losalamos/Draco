#!/usr/bin/env python

# Python module for interfacing with the module command on ASC
# systems.

import os, Utils

module_paths = os.environ['MODULEPATH'].split(':')
module_home  = os.environ['MODULESHOME']

def module_command(command):

    """module_command will execute the given module command"""

    exec os.popen('/usr/bin/modulecmd python %s' % command).read()
    

def get_avail_modules():

    """Return list of the available modules. Look in standard
    locations for recognized platforms. Look os.environ['MODULEPATH']

    flash: /usr/share/modules/modulefiles and subdirectories
    qsc:   /opt/modulefiles
    
    """

    modules = []

    paths = [path for path in module_paths if
             path.count('modulefiles') > 0]

    for path in paths:
        modules.append(Utils.listFiles(path, recurse=1))

        
    return modules

    
def get_loaded_modules():
    """Get the list of modules currently loaded."""

    return os.environ['LOADEDMODULES'].split(':')

    

if __name__=='__main__':

    module_command('remove tools/idl-6.1')
    module_command('list')
    module_command('add tools/idl-6.1')
    print get_loaded_modules()
    print get_avail_modules()


