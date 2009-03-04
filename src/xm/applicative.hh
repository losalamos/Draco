//----------------------------------*-C++-*----------------------------------//
// Copyright 1996 The Regents of the University of California. 
// All rights reserved.
//---------------------------------------------------------------------------//

#ifndef __xm_applicative_hh__
#define __xm_applicative_hh__

#ifndef __xm_xm_hh__
#error "Users should only include xm/xm.hh"
#endif

#include <cmath>
#include <cstdlib>

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
// Mathematical unary operators +,-.
//---------------------------------------------------------------------------//

//===========================================================================//
// class OpUnaryAdd - 

// 
//===========================================================================//

template<class P>
class OpUnaryAdd {
  public:
    static inline P apply( const P& x )
    {
	return +x;
    }
};

//===========================================================================//
// class OpUnarySub - 

// 
//===========================================================================//

template<class P>
class OpUnarySub {
  public:
    static inline P apply( const P& x )
    {
	return -x;
    }
};

//---------------------------------------------------------------------------//
// Mathematical intrinsic unary functions.
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
	return std::sin(x);
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
	return std::cos(x);
    }
};

//===========================================================================//
// class OpTan - 

// 
//===========================================================================//

template<class P>
class OpTan {
  public:
    static inline P apply( const P& x )
    {
	return std::tan(x);
    }
};

//===========================================================================//
// class OpACos - 

// 
//===========================================================================//

template<class P>
class OpACos {
  public:
    static inline P apply( const P& x )
    {
	return std::acos(x);
    }
};

//===========================================================================//
// class OpASin - 

// 
//===========================================================================//

template<class P>
class OpASin {
  public:
    static inline P apply( const P& x )
    {
	return std::asin(x);
    }
};

//===========================================================================//
// class OpATan - 

// 
//===========================================================================//

template<class P>
class OpATan {
  public:
    static inline P apply( const P& x )
    {
	return std::atan(x);
    }
};

//===========================================================================//
// class OpCosh - 

// 
//===========================================================================//

template<class P>
class OpCosh {
  public:
    static inline P apply( const P& x )
    {
	return std::cosh(x);
    }
};

//===========================================================================//
// class OpSinh - 

// 
//===========================================================================//

template<class P>
class OpSinh {
  public:
    static inline P apply( const P& x )
    {
	return std::sinh(x);
    }
};

//===========================================================================//
// class OpTanh - 

// 
//===========================================================================//

template<class P>
class OpTanh {
  public:
    static inline P apply( const P& x )
    {
	return std::tanh(x);
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
	return std::exp(x);
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
	return std::log(x);
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
	return std::log10(x);
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
	return std::sqrt(x);
    }
};

//===========================================================================//
// class OpCeil - 

// 
//===========================================================================//

template<class P>
class OpCeil {
  public:
    static inline P apply( const P& x )
    {
	return std::ceil(x);
    }
};

//===========================================================================//
// class OpAbs - 

// 
//===========================================================================//

template<class P>
class OpAbs {
  public:
    static inline P apply( const P& x )
    {
	return ( x > 0 ? x : -x );
    }
};

//===========================================================================//
// class OpLabs - 

// 
//===========================================================================//

template<class P>
class OpLabs {
  public:
    static inline P apply( const P& x )
    {
	return std::labs(x);
    }
};

//===========================================================================//
// class OpFabs - 

// 
//===========================================================================//

template<class P>
class OpFabs {
  public:
    static inline P apply( const P& x )
    {
	return std::fabs(x);
    }
};

//===========================================================================//
// class OpFloor - 

// 
//===========================================================================//

template<class P>
class OpFloor {
  public:
    static inline P apply( const P& x )
    {
	return std::floor(x);
    }
};

//---------------------------------------------------------------------------//
// Mathematical intrinsic binary functions.
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
      return std::pow(x,y);
    }
};

//===========================================================================//
// class OpATan2 - 

// 
//===========================================================================//

template<class P>
class OpATan2 {
  public:
    static inline P apply( const P& x, const P& y )
    {
      return std::atan2(x,y);
    }
};

//===========================================================================//
// class OpFmod - 

// 
//===========================================================================//

template<class P>
class OpFmod {
  public:
    static inline P apply( const P& x, const P& y )
    {
      return std::fmod(x,y);
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
