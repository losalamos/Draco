#! /usr/bin/env python
##---------------------------------------------------------------------------##
# Module: Finger
#
# Functions for talking to the LANL LDAP server via finger.
#
##---------------------------------------------------------------------------##

import os

def finger_ldap(keyword):

    """ Pass keyword to the finger daemon process as: 'finger
    <keyword>@lanl.gov to query the LDAP server
    
    >>> print finger_ldap("buksas")[0].split()[1:3]
    ['Buksas', 'Michael']
    """

    command = "finger %s@lanl.gov" % keyword
    stream = os.popen(command)

    output = stream.readlines()

    return output

##---------------------------------------------------------------------------##
class person:

    import re, sys

    znum    = "\d{6}"                 # six digits
    word    = "\w+"                   # one or more chars.
    initial = "\w{1}\."               # one char, followed by '.' 
    group   = "\w{1,5}-\d{1,2}"       # WWWWW-DD or fewer chars
    group2  = "\w+/\w+"               # WWW/WWW contractors etc...
    mail    = "\w{1}\d{3}"            # WDDD, exactly
    phone   = "\d{3}-\d{3}-\d{4}"     # DDD-DDD-DDDD
    email   = "\w+@(\w+\.)+(\w+)"     # blah@(here.)there
    space   = "\s+"                   # one or more whitespace chars.
    
    match_string =   "(?P<znum>%s)"          % znum    + space \
                   + "(?P<lname>%s)"         % word    + space \
                   + "(?P<fname>%s)"         % word    + space \
                   + "(?P<minit>%s)"         % initial + space \
                   + "(\((?P<nick>(%s))\))?" % word    + space \
                   + "(?P<group>%s|%s)"      % (group, group2) + space \
                   + "(?P<mail>%s)"          % mail    + space \
                   + "((?P<phone>%s))"       % phone   + space \
                   + "((?P<email>%s))"       % email

    matcher = re.compile(match_string)

    # Extract the keys from match_string
    keys = [key[2:-1] for key in re.findall("P<\w+>", match_string)]
    
    def __init__(self, finger_data):
        import sys
        
        """Match finger_data aginst the compiled regex and set local
        variables

        WARNING: This operation silently dies on incomplete LADP
        records. This needs to be fixed.

        """

        for key in person.keys:
            setattr(self, key, None)

        if finger_data:
            match = person.matcher.match(finger_data)

        if match:
            self.__dict__.update(match.groupdict())
        else:
            sys.stderr.write("Could not parse record: %s" % finger_data)

    def disp(self, out = sys.stdout):
        out.write("%s %s %s %s " % (self.znum,  self.lname,
                                    self.fname, self.minit))

        if self.nick: out.write("(%s) " % self.nick)

        out.write("%s %s %s %s\n" % (self.group, self.mail,
                                     self.phone, self.email))
                  

    def pretty_print(self, out = sys.stdout):
        out.write("       Name:  %s %s %s\n" % (self.fname,
                                                self.minit,
                                                self.lname)) 
        if self.nick:
            out.write("   Nickname:  %s\n" % self.nick)
        out.write("   Z number:  %s\n" % self.znum)
        out.write("      group:  %s\n" % self.group)
        out.write("  mail stop:  %s\n" % self.mail)
        out.write("      phone:  %s\n" % self.phone)
        out.write("      email:  %s\n" % self.email)
        


##---------------------------------------------------------------------------##
def finger(word):
    """function finger

    Return a list of people objects obtained from finger_ldap and the
    given word

    >>> me = finger("Buksas")[0]
    >>> me['lname']
    'Buksas'

    """

    matches = finger_ldap(word)

    people = [person(data) for data in matches]

    return people

##---------------------------------------------------------------------------##
def _test():
    import doctest, Finger
    doctest.testmod(Finger)


##---------------------------------------------------------------------------##
if __name__=="__main__":
    _test()
