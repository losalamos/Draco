//----------------------------------*-C++-*----------------------------------//
// Copyright 1996 The Regents of the University of California. 
// All rights reserved.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Created on: Sat Oct 26 12:54:12 1996
// Created by: Geoffrey Furnish
// Also maintained by:
//
//---------------------------------------------------------------------------//

extern "C" {
#include "Python.h"
}
#include "pythonrun.h"
#include "arrayobject.h"

static PyObject *ErrorObject;
#define Adfile_Error(msg) PyErr_SetString(ErrorObject, msg )

// Try to make an adfile a first class Python type.

#include "util/ADFile.hh"

class pyadfile {
    PyObject_HEAD
    ADFile *f;
  public:
    pyadfile( String name, int mode, int nrecs );
    ~pyadfile() { delete f; }

    PyObject *Get_keys();
    PyObject *Locate( ADKey& adk, int wrap );
    PyObject *len( int e );
    PyObject *binread( int n );

// Low level read operations.

    PyObject *read_int();
    PyObject *read_float();
    PyObject *read_double();
    PyObject *read_doubles( int n );
};

static void pyadfile_dealloc( pyadfile *adf );
static PyObject *pyadfile_getattr( pyadfile *adf, char *name );

/* 
   DEFINE THE TYPE OBJECT FORWARD REFERENCE:

   Every Python type must be associated with a (single shared) type object.
   Since we need the methods of the objects to define the type object,
   we will not define it here, instead we will create a forward reference
   since we will need to refer to the type object real soon.
*/

// staticforward PyTypeObject pyadfiletype;
// staticforward/statichere stuff is lunacy.  Do it RIGHT!!!!

statichere PyTypeObject pyadfiletype = {
    PyObject_HEAD_INIT(&PyType_Type)
    0,
    "pyadfile",
//     sizeof(pyadfileobject),
    sizeof(pyadfile),
    0,
    (destructor)pyadfile_dealloc,     /*tp_dealloc*/
    0,                               /*tp_print*/
    (getattrfunc)pyadfile_getattr,    /*tp_getattr*/
    0,                               /*tp_setattr*/
    0,                               /*tp_compare*/
//  (reprfunc)bstream_repr,          /*tp_repr*/
    0, //(reprfunc)bstream_repr,          /*tp_repr*/
    0,                               /*tp_as_number*/
//     &bstream_as_sequence,            /*tp_as_sequence*/
    0,//&bstream_as_sequence,            /*tp_as_sequence*/
    0,                               /*tp_as_mapping*/
    0,                               /*tp_hash*/
};


/* (staticforward is a macro that defines a forward reference
    convention that works with lots of compilers and C variants). */

/* 
     DEFINE THE TYPE TEST MACRO:

     The bstreamobject type test macro:
     it is usually useful to define a type test macro for a python
     type. Here we define is_bstreamobject 
*/

#define is_pyadfileobject(op) ((op)->ob_type == &pyadfiletype)

/* 
   INTERNAL UTILITY ROUTINES:

   Functions in this section are only used within this module.
   They are not exposed to Python or elsewhere.
*/

/***** METHODS ******/

pyadfile::pyadfile( String name, int mode, int nrecs )
{
    ob_type = &pyadfiletype;

    f = new ADFile( name, mode, nrecs );

    _Py_NewReference( this );
}

PyObject *pyadfile::Get_keys()
{
    Mat1<ADKey> keys = f->Get_keys();
    int nkeys = keys.size();

    PyObject *l = PyList_New(nkeys);
    for( int i=0; i < nkeys; i++ ) {
	ADKey key = keys[i];
	PyList_SetItem( l, i, Py_BuildValue("s", key.s) );
    }

    return l;
}

PyObject *pyadfile::Locate( ADKey& adk, int wrap )
{
    int e = f->Locate( adk, wrap );
    return Py_BuildValue("i", e );
}

PyObject *pyadfile::len( int e )
{
    return Py_BuildValue("i", f->len(e) );
}

PyObject *pyadfile::binread( int n )
{
    char *buf = new char[ n ];

    f->read( buf, n );

    PyObject *result = Py_BuildValue("s#", buf, n );
    delete[] buf;
    return result;
}

PyObject *pyadfile::read_int()
{
    int x;
    f->read( &x, sizeof(int) );
    return Py_BuildValue("i", x );
}

PyObject *pyadfile::read_float()
{
    float x;
    f->read( &x, sizeof(float) );
    return Py_BuildValue("f", x );
}

PyObject *pyadfile::read_double()
{
    double x;
    f->read( &x, sizeof(double) );
    return Py_BuildValue("d", x );
}

PyObject *pyadfile::read_doubles( int n )
{
    PyObject *ap = PyArray_FromDims( 1, &n, PyArray_DOUBLE );

    f->read( ((PyArrayObject *) ap)->data, n * sizeof(double) );

    return ap;
}

/*** STANDARD METHODS OF BSTREAM OBJECTS ***/

// Dealocate a pyadfile object.

static
void pyadfile_dealloc( pyadfile *adf )
{
    delete adf;
}

/**** NAMED (ATTRIBUTE) METHODS OF BSTREAM OBJECTS ****/

static
PyObject *pyadfile_Get_keys( pyadfile *self, PyObject *args )
{
    return self->Get_keys();
}

static
PyObject *pyadfile_Locate( pyadfile *self, PyObject *args )
{
    ADKey adk;
    char *key;
    int wrap;
    int rc = PyArg_ParseTuple( args, "si", &key, &wrap );
    if (!rc) return NULL;
    sprintf( adk.s, "%s", key );
    return self->Locate( adk, wrap );
}

static
PyObject *pyadfile_len( pyadfile *self, PyObject *args )
{
    int e;
    int rc = PyArg_ParseTuple( args, "i", &e );
    if (!rc) return NULL;
    return self->len(e);
}

static
PyObject *pyadfile_binread( pyadfile *self, PyObject *args )
{
    int n;
    int rc = PyArg_ParseTuple( args, "i", &n );
    if (!rc) return NULL;
    return self->binread(n);
}

static
PyObject *pyadfile_read_int( pyadfile *self, PyObject *args )
{
    return self->read_int();
}

static
PyObject *pyadfile_read_float( pyadfile *self, PyObject *args )
{
    return self->read_float();
}

static
PyObject *pyadfile_read_double( pyadfile *self, PyObject *args )
{
    return self->read_double();
}

static
PyObject *pyadfile_read_doubles( pyadfile *self, PyObject *args )
{
    int n;
    int rc = PyArg_ParseTuple( args, "i", &n );
    if (!rc) return NULL;
    return self->read_doubles( n );
}

/**** SEQUENCE METHODS ****/

// No sequencing support for adfile's.

static PyMethodDef pyadfile_methods[] = {

/* python name         (C type) C name               flag (always 1) */

    {"Get_keys",	(PyCFunction) pyadfile_Get_keys,	1},
    {"Locate",		(PyCFunction) pyadfile_Locate,		1},
    {"len",		(PyCFunction) pyadfile_len,		1},
    {"binread",		(PyCFunction) pyadfile_binread,		1},
    {"read_int",	(PyCFunction) pyadfile_read_int,	1},
    {"read_float",	(PyCFunction) pyadfile_read_float,	1},
    {"read_double",	(PyCFunction) pyadfile_read_double,	1},
    {"read_doubles",	(PyCFunction) pyadfile_read_doubles,	1},

  {NULL,           NULL}        /* sentinel */
};

/* THE STANDARD GETATTR FOR NAMED METHODS */

static
PyObject *pyadfile_getattr( pyadfile *adf, char *name )
{
    return Py_FindMethod( pyadfile_methods, (PyObject *) adf, name );
}


/* THE SEQUENCE METHODS STRUCTURE FOR ALL BSTREAM OBJECTS */

// static PySequenceMethods bstream_as_sequence = {
//   (inquiry)bstream_len,            /*sq_length*/
//   (binaryfunc)bstream_concat,      /*sq_concat*/
//   (intargfunc)bstream_repeat,      /*sq_repeat*/
//   (intargfunc)bstream_item,        /*sq_item*/
//   (intintargfunc)bstream_slice,    /*sq_slice*/
//   0,                               /*sq_ass_item*/
//   0,                               /*sq_ass_slice*/
// };

/* THE TYPEOBJECT FOR ALL BSTREAM OBJECTS (REFERS TO SEQUENCE METHODS) */

// statichere PyTypeObject bstreamtype = {
//   PyObject_HEAD_INIT(&PyType_Type)
//   0,
//   "bstream",
//   sizeof(bstreamobject),
//   0,
//   (destructor)bstream_dealloc,     /*tp_dealloc*/
//   0,                               /*tp_print*/
//   (getattrfunc)bstream_getattr,    /*tp_getattr*/
//   0,                               /*tp_setattr*/
//   0,                               /*tp_compare*/
//   (reprfunc)bstream_repr,          /*tp_repr*/
//   0,                               /*tp_as_number*/
//   &bstream_as_sequence,            /*tp_as_sequence*/
//   0,                               /*tp_as_mapping*/
//   0,                               /*tp_hash*/
// };

/***** MODULE METHODS *****/

static char *new_adfile_doc = "routine to make a python adfile";

extern "C"
PyObject *new_adfile( PyObject *self, PyObject *args )
{
// Parse args.
    char *name;
    int mode, nrecs;

    int rc = PyArg_ParseTuple( args, "sii", &name, &mode, &nrecs );
    if (!rc) { return NULL; }	// botched the call.

    pyadfile *pdf;

    try {
	pdf = new pyadfile( name, mode, nrecs );
    } catch(...) {
	Adfile_Error( "unable to make file." );

	Py_XINCREF(Py_None);
	return Py_None;
    }

    return (PyObject *) pdf;
}

/* MODULE METHODS STRUCTURE */

static struct PyMethodDef adfile_methods[] = {
    {"new_adfile", new_adfile, 1, new_adfile_doc},

    {NULL, NULL, NULL, NULL}		/* sentinel */
};

// Initialization function for the module, /must/ be called initadfile().

static char adfile_module_documentation[] = 
"Provides Python access to the ADFile utility.";

extern "C"
void initadfile()
{
// To set default value of global and static objects.
// Create the module and add the functions.

    cout << "executing initadfile" << endl;

    PyObject *m = Py_InitModule4( "adfile", adfile_methods,
				  adfile_module_documentation,
				  (PyObject *) NULL, PYTHON_API_VERSION );

// Add some symbolic constants to the module

    PyObject *d = PyModule_GetDict(m);
    ErrorObject = PyString_FromString( "adfile.error" );
    PyDict_SetItemString( d, "error", ErrorObject );

// XXXX Add constants here.

// Check for errors.
    if ( PyErr_Occurred() )
	Py_FatalError( "can't initialize module adfile" );
}

//---------------------------------------------------------------------------//
//                              end of adfilemodule.cc
//---------------------------------------------------------------------------//
