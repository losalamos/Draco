//----------------------------------*-C++-*----------------------------------//
// Copyright 1996 The Regents of the University of California. 
// All rights reserved.
//---------------------------------------------------------------------------//

#ifndef Py_Module_hh
#define Py_Module_hh

#include "Python.h"

#include "PyX_Objects.hh"

namespace Py {

//===========================================================================//
// class Module

// The purpose of this class is to serve as the C++ replacement for the
//===========================================================================//

    class Module : public PyObject {
	Module();
      public:
	friend Module *new_Module( char *name );

	Dict GetDict() {
	    PyObject *p = PyModule_GetDict (this);
	    Py_INCREF (p);
	    return Dict(p); 
	}
    };

    Module *new_Module( char *name );
}

Py::Module *XInitmodule4( char *name, struct PyMethodDef *methods, char *doc,
			  PyObject *passthrough, int module_api_version );

#endif				// Py_Module_hh

//---------------------------------------------------------------------------//
//                              end of Module.hh
//---------------------------------------------------------------------------//
