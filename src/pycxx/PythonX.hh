//----------------------------------*-C++-*----------------------------------//
// Copyright 1996 The Regents of the University of California. 
// All rights reserved.
//---------------------------------------------------------------------------//

#ifndef Py_PythonX_hh
#define Py_PythonX_hh

// Pick up the C API.
#include "Python.h"

// Now get the C++ PyArg class.
#include "PyX_Args.hh"
// Python objects in C++
#include "PyX_Objects.hh"
// Module support
#include "PyX_Module.hh"
#include <list>
// using std::list;

namespace Py {

    struct cxxmethodobject : public PyObject {
	struct PyMethodDef *m_ml;
	PyObject *m_self;
    };

    template<class T> struct X_MethodDef
    {
	typedef PyObject *(T::*mptr)( PyArgs *args );
	char *ml_name;
	mptr ml_meth;
	int ml_flags;
	char *ml_doc;
    };

    struct cxxmemberfuncobject : public PyObject {
	struct X_MethodDef<PyObject> *m_ml;
	PyObject *m_self;
    };
    extern PyTypeObject CxxMemberfunctype;

    PyObject *CxxMemberFunction_New( X_MethodDef<PyObject> *ml,
				     PyObject *self );

    PyObject *cxx_trampoline_mbrfnc_va( cxxmemberfuncobject *self, PyObject *args );
    PyObject *cxx_trampoline_mbrfnc_kw( cxxmemberfuncobject *self, PyObject *args );

    extern PyMethodDef tramp_mbrfnc_va;
    extern PyMethodDef tramp_mbrfnc_kw;

    template<class T> struct X_methodchain {
	X_MethodDef<T> *methods;
	X_methodchain *link;
    };

    template<class T>
	PyObject *CxxMemberFunction_New( X_MethodDef<T> *ml, T *self ) {

	cxxmemberfuncobject *op = PyObject_NEW( cxxmemberfuncobject,
						&CxxMemberfunctype );

	if (op) {
	// We have to have our own copy, since not type compatible.  The type
	// deallocator frees it when this object is collected.

	    X_MethodDef<PyObject> *x = (X_MethodDef<PyObject> *)
		malloc( sizeof(X_MethodDef<PyObject>) );

	    x->ml_name = ml->ml_name;
	    x->ml_meth = static_cast< X_MethodDef<PyObject>::mptr >
		( ml->ml_meth );
	    x->ml_flags = ml->ml_flags;
	    x->ml_doc = ml->ml_doc;	
    
	    op->m_ml = x;

	    Py_XINCREF(self);
	    op->m_self = self;
	}
	else
	    PyErr_SetString( PyExc_RuntimeError,
			     "Failed to allocate cxxmemberfuncobject." );

	return op;
    }

    template<class T>
	PyObject *newmemberfuncobject( X_MethodDef<T> *ml, T *self )
	{
	    PyObject *c = CxxMemberFunction_New( ml, self );
	    PyObject *v = PyCFunction_New( &tramp_mbrfnc_va, c );	
	    if (!c || !v) return NULL;
	    return v;
	}
    
    template<class T>
	PyObject *findmemberfuncinchain( X_methodchain<T> *chain,
					 T *self, char *name )
	{
	    if (strcmp(name, "__methods__") == 0) {
		PyErr_SetString( PyExc_RuntimeError,
				 "Too tired for this one, Bud." );
		return NULL;
	    // return listmethodchain(chain);
	    }
	    while (chain != NULL) {
		X_MethodDef<T> *ml = chain->methods;
		for (; ml->ml_name != NULL; ml++) {
		    if (name[0] == ml->ml_name[0] &&
			strcmp(name+1, ml->ml_name+1) == 0) {
			if (PyErr_Occurred()) {
			// We just /assume/ this was the failed attribute
			// lookup from a preceding base class.  This should
			// be checked in greater detail...
			    PyErr_Clear();
			}
			return newmemberfuncobject( ml, self );
		    }
		}
		chain = chain->link;
	    }
	    PyErr_SetString( PyExc_AttributeError, name );
	    return NULL;
	}

    typedef std::list<getattrfunc> parentlist;

/* Find a method in a single method list */

    template<class T>
	PyObject *FindMemberFunc( X_MethodDef<T> *methods, T *self,
				  char *name )
	{
	    struct X_methodchain<T> chain;
	    chain.methods = methods;
	    chain.link = NULL;
	    return findmemberfuncinchain( &chain, self, name );
	}

/* Search for a method first in a methodlist for a given class, and if not
   found, search the base classes for the method. */

    template<class T>
	PyObject *FindMemberFunc( X_MethodDef<T> *methods, T *self,
				  char *name, const parentlist& pl )
	{
	    struct X_methodchain<T> chain;
	    chain.methods = methods;
	    chain.link = NULL;
	    PyObject *result = findmemberfuncinchain( &chain, self, name );
	    if ( result == NULL )
		for( parentlist::const_iterator p = pl.begin();
		     p != pl.end(); p++ )
		    if ( result = (*p)( self, name ) ) break;
	    return result;
	}

    class CxxTypeBase : public PyObject {
      public:
	typedef PyObject *(PyObject::*getattrmeth)( char *name );

      protected:
	PyObject *pyo;
      public:
	CxxTypeBase( PyObject *_pyo );
	virtual ~CxxTypeBase() { cout << "deleting CxxTypeBase object" <<
				     endl; }
	virtual PyObject *_getattr( char *name ) =0;
    };

    template<class T> class CxxType
	: public CxxTypeBase
    {
	T *p;
	getattrfunc g;
	std::list<getattrfunc> parents;

      public:
	CxxType( T *_p, getattrfunc _g )
	    : CxxTypeBase( _p ), p(_p), g(_g) {}
	
	~CxxType() { delete p; }

	PyObject *_getattr( char *name )
	{
	    PyObject *res = (*g)( pyo, name );
	    return res;
	}
    };

}

#define Instantiate_helpers(T) \
template PyObject *Py::FindMemberFunc( Py::X_MethodDef< T > *, T *, char * ); \
template PyObject *Py::FindMemberFunc( Py::X_MethodDef< T > *, T *, char *, const parentlist& ); \
template PyObject *Py::findmemberfuncinchain( Py::X_methodchain< T > *, T *, char *); \
template PyObject *Py::newmemberfuncobject( Py::X_MethodDef< T > *ml, T *self ); \
template PyObject *Py::CxxMemberFunction_New( X_MethodDef< T > *ml, T *self );

#endif				// Py_PythonX_hh

//---------------------------------------------------------------------------//
//                              end of PythonX.hh
//---------------------------------------------------------------------------//
