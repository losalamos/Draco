//----------------------------------*-C++-*----------------------------------//
// XYZ_Quadratic.hh
// Geoffrey M. Furnish
// Wed Nov 19 17:13:40 1997
//---------------------------------------------------------------------------//
// @> 3d cartesian test case.
//---------------------------------------------------------------------------//

#ifndef __3T_test_XYZ_Quadratic_hh__
#define __3T_test_XYZ_Quadratic_hh__

#include "Quad_Params.hh"

//===========================================================================//
// class XYZ_Quadratic - Implement 3-d test case using quadratic forms

// This class implements a 3-d test case for the 3T package, using xyz
// cartesian geometry, and quadratic functional forms which can be correctly
// differenced to machine accuracy.
//===========================================================================//

class XYZ_Quadratic : private Quad_Params
{

  public:
    typedef Quad_Params params;

    XYZ_Quadratic( const Quad_Params& q ) : Quad_Params(q) {}

    double Et( double t ) const { return et0 + et1*t; }
    double Ex( double x ) const { return ex0 + ex1*x + ex2*x*x; }
    double Ey( double y ) const { return ey0 + ey1*y + ey2*y*y; }
    double Ez( double z ) const { return ez0 + ez1*z + ez2*z*z; }

    double Exyz( double x, double y, double z ) const
    {
	return Ex(x)*Ey(y)*Ez(z);
    }
    double E( double x, double y, double z, double t ) const
    {
	return Et(t) * Ex(x) * Ey(y) * Ez(z);
    }

    double dEtdt() const { return et1; }
    double dExdx( double x ) const { return ex1 + 2. * ex2 * x; }
    double dEydy( double y ) const { return ey1 + 2. * ey2 * y; }
    double dEzdz( double z ) const { return ez1 + 2. * ez2 * z; }

    double d2Exdxx( double x ) const { return 2. * ex2; }
    double d2Eydyy( double y ) const { return 2. * ey2; }
    double d2Ezdzz( double z ) const { return 2. * ez2; }

    double Dx( double x ) const { return dx0 + dx1*x; }
    double Dy( double y ) const { return dy0 + dy1*y; }
    double Dz( double z ) const { return dz0 + dz1*z; }

    double D( double x, double y, double z ) const
    {
	return Dx(x) * Dy(y) * Dz(z);
    }

    double dDxdx() const { return dx1; }
    double dDydy() const { return dy1; }
    double dDzdz() const { return dz1; }
};

#endif                          // __3T_test_XYZ_Quadratic_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/test/XYZ_Quadratic.hh
//---------------------------------------------------------------------------//
