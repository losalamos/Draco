import array

# Functions for mappings:
def add(a,b): return a+b
def mult(a,b): return a*b

def AddArrays(a,b):
    "Return the element-wise addition of two arrays."
    c = map(add,a,b)
    return array.array(a.typecode,c)

def Accumulate(a,b):
    "Add the elements in the second argument to the first"
    for i in range(len(b)):
        a[i] += b[i]

def AccumulateSquare(a,b):
    "Add the square of the elements in the second argument to the first."
    for i in range(len(b)):
        a[i] += b[i]**2

def AccumulateFunc(a, b, func):
    "Add func applied to the elements of the second argument to the first."
    for i in range(len(b)):
        a[i] += func(b[i])
    
def Scale(a,beta):
    "Multiple the elements of the first argument by the scalar second argument"
    for i in range(len(a)):
        a[i] *= beta

def ComputeSum(arg):
    "Return the sum of the elements in an array."
    reduce(add, arg)

def ComputeMean(arg):
    "Return the mean value of an array"
    return ComputeSum(arg)/len(arg)       