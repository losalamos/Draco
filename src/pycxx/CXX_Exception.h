//----------------------------------*-C++-*----------------------------------//
// Copyright 1998 The Regents of the University of California. 
// All rights reserved. See LEGAL.LLNL for full text and disclaimer.
//---------------------------------------------------------------------------//

#ifndef __CXX_Exception_h
#define __CXX_Exception_h

#include "Python.h"

#include <string>
#include <iostream>

namespace Py {
    
    class Exception {
    public:
        explicit Exception () {}
        
        Exception (const std::string& reason) {
            PyErr_SetString (PyExc_RuntimeError, reason.c_str());
        }
        
        Exception (PyObject* exception, const std::string& reason) {
            PyErr_SetString (exception, reason.c_str());
        }
        
        void clear() // clear the error
            // technically but not philosophically const
        {
            PyErr_Clear();
        }
    };
    
    class TypeError: public Exception  {
    public:
        TypeError (const std::string& reason)
            : Exception()
        {
            PyErr_SetString (PyExc_TypeError,reason.c_str());
        }
    };
    
    class IndexError: public Exception  {
    public:
        IndexError (const std::string& reason)
            : Exception()
        {
            PyErr_SetString (PyExc_IndexError, reason.c_str());
        }
    };
    
    class AttributeError: public Exception {
    public:
        AttributeError (const std::string& reason)
            : Exception()
        {
            PyErr_SetString (PyExc_AttributeError, reason.c_str());
        }		
    };
    
    class NameError: public Exception  {
    public:
        NameError (const std::string& reason)
            : Exception()
        {
            PyErr_SetString (PyExc_NameError, reason.c_str());
        }
    };
    
    class RuntimeError: public Exception	{
    public:
        RuntimeError (const std::string& reason)
            : Exception()
        {
            PyErr_SetString (PyExc_RuntimeError, reason.c_str());
        }
    };
    
    class SystemError: public Exception  {
    public:
        SystemError (const std::string& reason)
            : Exception()
        {
            PyErr_SetString (PyExc_SystemError,reason.c_str());
        }
    };
    
    class KeyError: public Exception	{
    public:
        KeyError (const std::string& reason)
            : Exception()
        {
            PyErr_SetString (PyExc_KeyError,reason.c_str());
        }
    };
    
    
    class ValueError: public Exception {
    public:
        ValueError (const std::string& reason)
            : Exception()
        {
            PyErr_SetString (PyExc_ValueError, reason.c_str());
        }
    };
    
    class OverflowError: public Exception  {
    public:
        OverflowError (const std::string& reason)
            : Exception()
        {
            PyErr_SetString (PyExc_OverflowError, reason.c_str());
        }		
    };
    
    class ZeroDivisionError: public Exception  {
    public:
        ZeroDivisionError (const std::string& reason)
            : Exception() 
        {
            PyErr_SetString (PyExc_ZeroDivisionError, reason.c_str());
        }
    };
    
    class MemoryError: public Exception  {
    public:
        MemoryError (const std::string& reason)
            : Exception()
        {
            PyErr_SetString (PyExc_MemoryError, reason.c_str());
        }	
    };
    
    class SystemExit: public Exception  {
    public:
        SystemExit (const std::string& reason)
            : Exception() 
        {
            PyErr_SetString (PyExc_SystemExit,reason.c_str());
        }
    };
    
} // end namespace Py

#endif
