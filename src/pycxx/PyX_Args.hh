//----------------------------------*-C++-*----------------------------------//
// Copyright 1996 The Regents of the University of California. 
// All rights reserved.
//---------------------------------------------------------------------------//

#ifndef PyX_Args_hh
#define PyX_Args_hh

#include "Python.h"

//#include "PyX_xcpt.hh"
#include "PyX_Objects.hh"

#include <string>

namespace Py {

template<class T>
class pytraits {
  public:
    typedef PyObject *holder_t;
    static const char *fmt() { return "O"; }
    static void post_config( const holder_t& h, T& t ) { t = h; }
};

template<> class pytraits<int> {
  public:
    typedef int holder_t;
    static const char *fmt() { return "i"; }
    static void post_config( const holder_t& h, int& t ) { t = h; }
};

template<> class pytraits<float> {
  public:
    typedef float holder_t;
    static const char *fmt() { return "f"; }
    static void post_config( const holder_t& h, float& t ) { t = h; }
};

template<> class pytraits<double> {
  public:
    typedef double holder_t;
    static const char *fmt() { return "d"; }
    static void post_config( const holder_t& h, double& t ) { t = h; }
};

//===========================================================================//
// class PyArgs

// The purpose of this class is to serve as the C++ replacement for the
// "PyObject *args" which is the conventional second argument to compiled
// extension functions. 
//===========================================================================//

class PyArgs : public PyObject {
  public:
    int ParseTuple()
    {
	int result = PyArg_ParseTuple( this, "" );
	if ( !result ) {
	    throw PyException();
	}
	return result;
    }

    template<class X1>
    int ParseTuple( X1& x1 )
    {
	pytraits<X1>::holder_t h1;

	std::string fmt = pytraits<X1>::fmt();

	int result = PyArg_ParseTuple( this,
				       const_cast<char *>(fmt.data()),
				       &h1 );
	if (!result) throw PyException();

	pytraits<X1>::post_config( h1, x1 );

	return result;
    }

    template<class X1, class X2>
    int ParseTuple( X1& x1, X2& x2 )
    {
	pytraits<X1>::holder_t h1;
	pytraits<X2>::holder_t h2;

	std::string fmt = pytraits<X1>::fmt();
	fmt += pytraits<X2>::fmt();

	int result = PyArg_ParseTuple( this,
				       const_cast<char *>(fmt.data()),
				       &h1, &h2 );
	if (!result) throw PyException();

	pytraits<X1>::post_config( h1, x1 );
	pytraits<X2>::post_config( h2, x2 );

	return result;
    }

    template<class X1, class X2, class X3>
    int ParseTuple( X1& x1, X2& x2, X3& x3 )
    {
	pytraits<X1>::holder_t h1;
	pytraits<X2>::holder_t h2;
	pytraits<X3>::holder_t h3;

	std::string fmt = pytraits<X1>::fmt();
	fmt += pytraits<X2>::fmt();
	fmt += pytraits<X3>::fmt();

	int result = PyArg_ParseTuple( this,
				       const_cast<char *>(fmt.data()),
				       &h1, &h2, &h3 );
	if (!result) throw PyException();

	pytraits<X1>::post_config( h1, x1 );
	pytraits<X2>::post_config( h2, x2 );
	pytraits<X3>::post_config( h3, x3 );

	return result;
    }
};

}

#endif				// PyX_Args_hh

//---------------------------------------------------------------------------//
//                              end of PyX_Args.hh
//---------------------------------------------------------------------------//
