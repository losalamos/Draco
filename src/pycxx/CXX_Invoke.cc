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
#if 0
/* Methods (the standard built-in methods, that is) */

    template<class FS>
    void cxx_meth_dealloc( cxxmethodobject<FS> *m )
    {
        Py_XDECREF( m->self );
        free( (char *) m );
    }

    template<class FS>
    PyObject *cxx_meth_getattr( cxxmethodobject<FS> *m, char *name )
    {
        if (strcmp(name, "__name__") == 0) {
            return newstringobject(m->md->name);
        }
        if (strcmp(name, "__doc__") == 0) {
            char *doc = m->m_ml->ml_doc;
            if (doc != NULL)
                return newstringobject(doc);
            INCREF(None);
            return None;
        }
        if (strcmp(name, "__self__") == 0) {
            object *self;
            if (getrestricted()) {
                err_setstr(RuntimeError,
                           "method.__self__ not accessible in restricted mode");
                return NULL;
            }
            self = m->self;
            if (self == NULL)
                self = None;
            INCREF(self);
            return self;
        }
        if (strcmp(name, "__members__") == 0) {
            return mkvalue("[sss]", "__doc__", "__name__", "__self__");
        }
        err_setstr(AttributeError, name);
        return NULL;
    }

    template<class FS>
    PyObject *cxx_meth_repr( cxxmethodobject<FS> *m )
    {
        char buf[200];
        if (m->m_self == NULL)
            sprintf(buf, "<built-in function %.80s>", m->md->name);
        else
            sprintf(buf,
                    "<built-in method %.80s of %.80s object at %lx>",
                    m->md->name, m->self->ob_type->tp_name,
                    (long)m->m_self);
        return newstringobject(buf);
    }

    template<class FS>
    int cxx_meth_compare( cxxmethodobject<FS> *a, cxxmethodobject<FS> *b )
    {
        if (a->m_self != b->m_self)
            return cmpobject(a->m_self, b->m_self);
        if (a->m_ml->ml_meth == b->m_ml->ml_meth)
            return 0;
        if (strcmp(a->m_ml->ml_name, b->m_ml->ml_name) < 0)
            return -1;
        else
            return 1;
    }

    template<class FS>
    long cxx_meth_hash( cxxmethodobject<FS> *a )
    {
        long x;
        if (a->m_self == NULL)
            x = 0;
        else {
            x = hashobject(a->m_self);
            if (x == -1)
                return -1;
        }
        return x ^ (long) a->m_ml->ml_meth;
    }
#endif
#if 0
    typeobject CxxMethodtype = {
        OB_HEAD_INIT(&Typetype)
        0,
        "cxx_builtin_function_or_method",
        sizeof(cxxmethodobject),
        0,
        (destructor)cxx_meth_dealloc, /*tp_dealloc*/
        0,		/*tp_print*/
        (getattrfunc)cxx_meth_getattr, /*tp_getattr*/
        0,		/*tp_setattr*/
        (cmpfunc)cxx_meth_compare, /*tp_compare*/
        (reprfunc)cxx_meth_repr, /*tp_repr*/
        0,		/*tp_as_number*/
        0,		/*tp_as_sequence*/
        0,		/*tp_as_mapping*/
        (hashfunc)cxx_meth_hash, /*tp_hash*/
    };
#endif
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
        //            return (*self->m_md->meth)( self->m_self, args );
//         //            return trampoline_traits<FS>::invoke( self->self, args );8
//             throw Exception( "don't know how to invoke a method yet." );
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
//         cout << "Bad trouble!" << endl;
//         return Nothing();
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
//         Module *m;
        PyObject *m;
//         Dict d;
        PyObject *v;
//         PyMethodDef *ml;
        if (module_api_version != PYTHON_API_VERSION)
            fprintf( stderr, api_version_warning,
                     name, PYTHON_API_VERSION, name, module_api_version);
        try {
//             m = new_Module( name );
            m = PyImport_AddModule( name );
            if (!m) throw Exception();
        }
        catch( Exception &pyx ) {
            fprintf(stderr, "initializing module: %s\n", name);
            Py_FatalError("can't create a module");
        }
//         d = m->GetDict();
        Dict d( PyModule_GetDict(m) );

    // Subterfuge.  Okay, what we are going to actually do is register a
    // callback to the C++ module function dispatch routine, and pass it the
    // methoddef struct as _its_ clientdata, which will allow it to call the
    // module's C++ function.

        for( List::iterator method = methods.begin();
             method != methods.end();
             method++ )
        {
            cout << "trying to register a method." << endl;
//         for( ml = methods; ml->ml_name != NULL; ml++ ) {
//             PyObject *c = CxxFunction_New( ml, passthrough );
            PyObject *c = CxxFunction_New( *method, passthrough );
            if (!c) cout << "CXX_Function_New flunked.\n";

            PyObject *v =
//                 PyCFunction_New( (PyCFunction) &cxx_trampoline_va<CXXFunction>, c );
                PyCFunction_New( &tramp_vad, c );
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
            if (v == NULL || PyDict_SetItemString(d.lend (), ml->ml_name, v) != 0) {
                fprintf(stderr, "initializing module: %s\n", name);
                Py_FatalError("can't initialize module");
            }
            Py_DECREF(v);
#endif

//             if (v == NULL || PyDict_SetItemString(d.lend (), ml->ml_name, v) != 0) {
            if (v == NULL || PyDict_SetItemString(d.ptr (), (*method)->name, v) != 0) {
                fprintf(stderr, "initializing module: %s\n", name);
                Py_FatalError("can't initialize module");
            }
            Py_DECREF(v);


        }
        if (module_doc != NULL) {
            v = PyString_FromString(module_doc);
//             if (v == NULL || PyDict_SetItemString(d.lend (), "__doc__", v) != 0)
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
