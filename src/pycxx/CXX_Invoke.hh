//----------------------------------*-C++-*----------------------------------//
// CXX_Invoke.hh
// Geoffrey M. Furnish
// Thu Apr  2 16:35:41 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __pycxx_CXX_Invoke_hh__
#define __pycxx_CXX_Invoke_hh__

#include "CXX_Objects.h"
#include "CXX_Extensions.h"

#include <assert.h>

namespace Py 
{
    typedef PyObject *(*CXXFunction) (Object, Tuple);
    typedef PyObject *(*CXXFunctionWithKeywords)
        (Object, Tuple *, Dict *);

    template<class FunctionSignature>
    struct CXXMethodDef {
        char *name;
        FunctionSignature meth;
        char *doc;

        typedef FunctionSignature FS;

//         CXXMethodDef() { assert(0); }
        CXXMethodDef( char *n, FunctionSignature m, char *d )
            : name(n), meth(m), doc(d) {}
    };

    template<class FunctionSignature>
    struct cxxmethodobject
        :
//         public PyObject
        public PythonExtension< cxxmethodobject<FunctionSignature> >
    {
        CXXMethodDef<FunctionSignature> *md;
        PyObject *self;

        cxxmethodobject( CXXMethodDef<FunctionSignature> *md_,
                         PyObject *s )
            : PythonExtension< cxxmethodobject<FunctionSignature> >(),
              md(md_),
              self(s)
        {}

        PyObject *invoke( PyObject *args )
        {
//             cout << "Con't know how to invoke yet." << endl;
//             return Nothing();
            return (*(md->meth))( Object(self), Tuple(args) );
        }
    };

    template<class FS>
    PyObject *cxx_trampoline_va( cxxmethodobject<FS> *self, PyObject *args );

    template<class List>
    PyObject *XInitmodule4( char *name, List& methods,
                            char *module_doc,
                            PyObject *passthrough, int module_api_version );
}

#endif                          // __pycxx_CXX_Invoke_hh__

//---------------------------------------------------------------------------//
//                              end of pycxx/CXX_Invoke.hh
//---------------------------------------------------------------------------//
