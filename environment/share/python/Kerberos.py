"""Utilities for working with the Kerberos authentication system.
"""

import os

klist    = "/usr/local/j2re1.4.2_05/bin/klist"
k5init   = "/usr/local/bin/k5init -f"
kdestroy = "/usr/kerberos/bin/kdestroy"

def has_ticket():
    """Check to see if the executing user has a valid kerberos ticket.

    Uses the return value of klist to perform the test.

    XXX This doesn't always work. klist sometimes returns true even
    though the last ticket has expired.

    Could parse the output from klist in this case. It doesn't say
    "expired" but gives an expiration time.

    """

    command = os.popen(klist, 'r')

    ret = command.close()

    return (ret == None)


def get_ticket():
    """Run k5init -f to get a kerberos ticket for the current user."""

    ret = os.system(k5init)

    return ret


def destroy_ticket():
    """Run kdestroy to eliminate any existing Kerberos tickets"""

    return os.system(kdestroy)


    


    
