//----------------------------------*-C++-*----------------------------------//
// Copyright 1996 The Regents of the University of California. 
// All rights reserved.
//---------------------------------------------------------------------------//

#ifndef __xm_applicative_hh__
#define __xm_applicative_hh__

#ifndef __xm_xm_hh__
#error "Users should only include xm/xm.hh"
#endif

#ifdef __sgi
#include <math.h>
#else
#include <cmath>
#endif

XM_NAMESPACE_BEG

//---------------------------------------------------------------------------//
// Applicative templates useful to all ranks since they apply only to
// elements of the containers which form the expression.
//---------------------------------------------------------------------------//
// Mathematical binary operators +,-,*,/.
//---------------------------------------------------------------------------//

//===========================================================================//
// class OpAdd - 

// 
//===========================================================================//

template<class P>
class OpAdd {
  public:
    static inline P apply( const P& a, const P& b )
    {
	return a + b;
    }
};

//===========================================================================//
// class OpSub - 

// 
//===========================================================================//

template<class P>
class OpSub {
  public:
    static inline P apply( const P& a, const P& b )
    {
	return a - b;
    }
};

//===========================================================================//
// class OpMul - 

// 
//===========================================================================//

template<class P>
class OpMul {
  public:
    static inline P apply( const P& a, const P& b )
    {
	return a * b;
    }
};

//===========================================================================//
// class OpDiv - 

// 
//===========================================================================//

template<class P>
class OpDiv {
  public:
    static inline P apply( const P& a, const P& b )
    {
	return a / b;
    }
};

//---------------------------------------------------------------------------//
// Mathematical intrinsic functions.
//---------------------------------------------------------------------------//

//===========================================================================//
// class OpSin - 

// 
//===========================================================================//

template<class P>
class OpSin {
  public:
    static inline P apply( const P& x )
    {
	return ::sin(x);
    }
};

//===========================================================================//
// class OpCos - 

// 
//===========================================================================//

template<class P>
class OpCos {
  public:
    static inline P apply( const P& x )
    {
	return ::cos(x);
    }
};

//===========================================================================//
// class OpExp - 

// 
//===========================================================================//

template<class P>
class OpExp {
  public:
    static inline P apply( const P& x )
    {
	return ::exp(x);
    }
};

//===========================================================================//
// class OpLog - 

// 
//===========================================================================//

template<class P>
class OpLog {
  public:
    static inline P apply( const P& x )
    {
	return ::log(x);
    }
};

//===========================================================================//
// class OpLog10 - 

// 
//===========================================================================//

template<class P>
class OpLog10 {
  public:
    static inline P apply( const P& x )
    {
	return ::log10(x);
    }
};

//===========================================================================//
// class OpSqrt - 

// 
//===========================================================================//

template<class P>
class OpSqrt {
  public:
    static inline P apply( const P& x )
    {
	return ::sqrt(x);
    }
};

//---------------------------------------------------------------------------//
// A special case (a binary mathematical intrinsic)...
//---------------------------------------------------------------------------//

//===========================================================================//
// class OpPow - 

// 
//===========================================================================//

template<class P>
class OpPow {
  public:
    static inline P apply( const P& x, const P& y )
    {
	return ::pow(x,y);
    }
};

//---------------------------------------------------------------------------//
// Relationals
//---------------------------------------------------------------------------//

//===========================================================================//
// class OpMin - 

// 
//===========================================================================//

template<class P>
class OpMin {
  public:
    static inline P apply( const P& x, const P& y )
    {
	return x < y ? x : y;
    }
};

//===========================================================================//
// class OpMax - 

// 
//===========================================================================//

template<class P>
class OpMax {
  public:
    static inline P apply( const P& x, const P& y )
    {
	return x > y ? x : y;
    }
};

XM_NAMESPACE_END

#endif                          // __xm_applicative_hh__

//---------------------------------------------------------------------------//
//                              end of xm/applicative.hh
//---------------------------------------------------------------------------//
