//----------------------------------*-C++-*----------------------------------//
// CXX_Invoke.cc
// Geoffrey M. Furnish
// Thu Apr  2 16:35:41 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "CXX_Invoke.hh"
#include "CXX_Exception.h"

namespace Py {
    template<class FS>
    PyObject *CxxFunction_New( CXXMethodDef<FS> *md, PyObject *self )
    {
        Py_XINCREF(self);
        cxxmethodobject<FS> *m = new cxxmethodobject<FS>( md, self );
        return static_cast<PyObject *>(m);
    }

    template<class FS>
    PyObject *cxx_trampoline_va( cxxmethodobject<FS> *self, PyObject *args )
    {
        try {
            return self->invoke( args );
        }
        catch( Exception& pyx ) {
        // This happens when the Python C API has already set an exception
        // condition, which is indicated by a return value of NULL.
            return NULL;
        }
        catch( const char *msg ) {
            char buf[ 1024 ];
            sprintf( buf, "a C++ exception has occurred: %s", msg );
            PyErr_SetString( PyExc_RuntimeError, buf );
            return NULL;
        }
        catch(...) {
            PyErr_SetString( PyExc_RuntimeError,
                             "an unrecognized C++ exception has occurred." );
            return NULL;
        }
    }

    PyObject *tramp_va( PyObject *self, PyObject *args )
    {
        cxxmethodobject<CXXFunction> *pf =
            static_cast<cxxmethodobject<CXXFunction> *>( self );
        if (!pf) cout << "Bad trouble!" << endl;
        return cxx_trampoline_va<CXXFunction>( pf, args );
    }

    PyMethodDef tramp_vad = { "tramp_va", (PyCFunction) tramp_va, 1, "xxx" };

//     PyObject *cxx_trampoline_kw( PyObject *self, PyObject *args, PyObject *kw )
//     {
//         PyErr_SetString( PyExc_RuntimeError,
//                          "Not able to invoke C++ keyword arg function yet." );
//         return NULL;
//     }

/* initmodule4() parameters:
   - name is the module name
   - methods is the list of top-level functions
   - doc is the documentation string
   - passthrough is passed as self to functions defined in the module
   - api_version is the value of PYTHON_API_VERSION at the time the
   module was compiled
*/

    static char api_version_warning[] =
    "WARNING: Python C API version mismatch for module %s:\n\
This Python has API version %d, module %s has version %d.\n";

    template<class List>
    PyObject *XInitmodule4( char *name, List& methods,
                            char *module_doc,
                            PyObject *passthrough, int module_api_version )
    {
        PyObject *m;
        PyObject *v;
        if (module_api_version != PYTHON_API_VERSION)
            fprintf( stderr, api_version_warning,
                     name, PYTHON_API_VERSION, name, module_api_version);

        try {
            m = PyImport_AddModule( name );
            if (!m) throw Exception();
        }
        catch( Exception &pyx ) {
            fprintf(stderr, "initializing module: %s\n", name);
            Py_FatalError("can't create a module");
        }
        Dict d( PyModule_GetDict(m) );

    // Subterfuge.  Okay, what we are going to actually do is register a
    // callback to the C++ module function dispatch routine, and pass it the
    // methoddef struct as _its_ clientdata, which will allow it to call the
    // module's C++ function.

        for( List::iterator method = methods.begin();
             method != methods.end();
             method++ )
        {
        // Here we build an object which is able to actually invoke the
        // C++ function. 

            PyObject *c = CxxFunction_New( *method, passthrough );
            if (!c) cout << "CXX_Function_New flunked.\n";

        // Here we build a C method object which runs a function that knows
        // how to ask the C++ method object to invoke the function it
        // represents. 

            PyObject *v = PyCFunction_New( &tramp_vad, c );
            if (!v) cout << "Can't make C func to invoke tramp_va" << endl;

#if 0
            if (ml->ml_flags == 1) {
                v = PyCFunction_New( &tramp_va, c );
            }
            else if (ml->ml_flags == 2) {
                v = PyCFunction_New( &tramp_kw, c );
            }
            else {
                fprintf(stderr, "initializing module: %s\n", name);
                Py_FatalError("can't initialize module");
            }
#endif

            if ( v == NULL ||
                 PyDict_SetItemString(d.ptr (), (*method)->name, v) != 0 ) {
                fprintf(stderr, "initializing module: %s\n", name);
                Py_FatalError("can't initialize module");
            }
            Py_DECREF(v);


        }
        if (module_doc != NULL) {
            v = PyString_FromString(module_doc);
            if (v == NULL || PyDict_SetItemString( d.ptr(), "__doc__", v) != 0)
                Py_FatalError("can't add doc string");
            Py_DECREF(v);
        }
        return m;
    }
}

//---------------------------------------------------------------------------//
//                              end of CXX_Invoke.cc
//---------------------------------------------------------------------------//
