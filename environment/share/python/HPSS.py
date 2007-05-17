"""Pacakge defined facilities for working with HPSS and the Mercury
system."""

import exceptions

import Verbosity, Utils

##---------------------------------------------------------------------------##
def store(file, comment = '', hpss_path = "xfer", verbose = Verbosity.ignore()):

    """Store a single file in users's HPSS space. Default directory:
    xfer.

    Also, accept a verbosity function for logging.

    """

    import os

    (file_path, file_name) = os.path.split(file)
    hpss_file = "%s/%s" % (hpss_path, file_name)
    
    verbose("Storing file %s in HPSS as %s" % (file, hpss_file))

    if comment:  comment_arg = "--cmt \"%s\"" % (comment,)
    else:        comment_arg = ''


    directory_arg = "-d %s" % (hpss_path,)

    command = ["/usr/local/bin/psi", "store",
               comment_arg,
               directory_arg,
               file_name] 
    
    verbose("... with command: %s" % " ".join(command), 2)

    os.chdir(file_path)
    store_status = os.spawnv(os.P_WAIT, command[0], command)
        
    verbose(" ... return status: %d" % (store_status,), 2)

    if store_status != 0:
        raise exceptions.RuntimeError("Unable to store file %s in HPSS. Psi command " \
                                      "return code is %s" % (file, store_status,))
                           
    return hpss_file



##---------------------------------------------------------------------------##
def push(hpss_filename, tag = '', verbose  = Verbosity.ignore()):

    """Push a function in HPSS via the Mercury system. User provides a
    filename and an optional id tag. Also accepts an optional
    verbosity function"""

    import os

    if tag:
        file_id = "id=%s" % (tag,)
    else:
        file_id = ''

    command = ["/usr/local/bin/push", file_id, hpss_filename]
        
    verbose("Pushing file %s" % (hpss_filename,))
    verbose("...with command: %s" % (" ".join(command),), 2)

    push_status = os.spawnv(os.P_WAIT, command[0], command)

    if push_status != 0:
        raise exceptions.RunTimeError("Failed to push file: %s with error " +
                                      "condition: %s " % (hpss_filename, push_status))




