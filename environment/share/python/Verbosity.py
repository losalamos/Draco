import sys

class Verb(object):

    def stdout(x) : print x
    stdout = staticmethod(stdout)

    def __init__(self, level = 1, function = stdout):

        self.level    = level
        self.function = function

    def __call__(self, message, level = 1):

        if (level <= self.level):
            self.function(message)

def make(l = 1, function = Verb.stdout): return Verb(l, function)

def all(function = Verb.stdout): return make(sys.maxint, function)

def ignore(): return make(0)
    
