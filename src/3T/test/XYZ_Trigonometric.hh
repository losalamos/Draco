//----------------------------------*-C++-*----------------------------------//
// XYZ_Trigonometric.hh
// Scott A. Turner (based on Geoffrey Furnish's XYZ_Quadratic.hh)
// 11 December 1997
//---------------------------------------------------------------------------//
// @> 3d cartesian test case.
//---------------------------------------------------------------------------//

#ifndef __3T_test_XYZ_Trigonometric_hh__
#define __3T_test_XYZ_Trigonometric_hh__

#include "Trig_Params.hh"

//===========================================================================//
// class XYZ_Trigonometric - Implement 3-d test case using trigonometric forms

// This class implements a 3-d test case for the 3T package, using xyz
// cartesian geometry, and trigonometric functional forms.
//===========================================================================//

class XYZ_Trigonometric : private Trig_Params
{

  public:
    typedef Trig_Params params;

    XYZ_Trigonometric( const Trig_Params& q ) : Trig_Params(q) {}

    double Et( double t ) const { return et0 + et1*t; }
    double Ex( double x ) const { return ex0 + ex1*sin(x); }
    double Ey( double y ) const { return ey0 + ey1*cos(y); }
    double Ez( double z ) const { return ez0 + ez1*tan(z); }

    double Exyz( double x, double y, double z ) const
    {
	return Ex(x) * Ey(y) * Ez(z);
    }
    double E( double x, double y, double z, double t ) const
    {
	return Et(t) * Ex(x) * Ey(y) * Ez(z);
    }

    double dEtdt() const { return et1; }
    double dExdx( double x ) const { return  ex1*cos(x); }
    double dEydy( double y ) const { return -ey1*sin(y); }
    double dEzdz( double z ) const { return  ez1/cos(z)/cos(z); }

    double d2Exdxx( double x ) const { return -ex1*sin(x); }
    double d2Eydyy( double y ) const { return -ey1*cos(y); }
    double d2Ezdzz( double z ) const { return  ez1*2.0/cos(z)/cos(z)*tan(z); }

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

#endif                          // __3T_test_XYZ_Trigonometric_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/test/XYZ_Trigonometric.hh
//---------------------------------------------------------------------------//
