//----------------------------------*-C++-*----------------------------------//
// Test_Prob.hh
// Geoffrey Furnish
// Tue Sep 30 16:52:41 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __t3_Test_Prob_hh__
#define __t3_Test_Prob_hh__

#include "t3/Run_DB.hh"

#include "c4/NodeInfo.hh"

#include "ds++/Mat.hh"

#include "linalg/pcg_DB.hh"

class ADFile;

//===========================================================================//
// class Test_Prob - 

// 
//===========================================================================//

class Test_Prob : private Run_DB, private C4::NodeInfo {
  public:

    int nct;			// # of total cells in problem.
    int ncp;			// # of cells on this processor.

    int goff;			// global offset of this processor's cells.
    int nxy;

// Convert a local cell index to its (i,j,k) indexes in whole domain.

    int I(int x) const { x += goff; return (x % nxy) % ncx; }
    int J(int x) const { x += goff; return (x % nxy) / ncx; }
    int K(int x) const { x += goff; return  x / nxy; }

    int goffset( int i, int j, int k ) const { return ncx*(k*ncy + j) +i; }
    int goffset( int nc ) const
    {
	nc += goff;
	return goffset( I(nc), J(nc), K(nc) );
    }

    Mat2<double> A;
    Mat1<double> xc, yc, zc;
    Mat1<double> xf, yf, zf;
    double       dx, dy, dz;
    Mat1<double> vc;
    Mat1<double> xA, yA, zA;

    pcg_DB       pcg_db;

    ADFile *adf;

  public:
    Test_Prob();
//     Test_Prob( const Test_Prob& );
    ~Test_Prob();
//     Test_Prob& operator=( const Test_Prob& );

    void run();

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
	return //(dx0 + dx1*x) * (dy0 + dy1*y) * (dz0 + dz1*z);
	    Dx(x) * Dy(y) * Dz(z);
    }

    double dDxdx() const { return dx1; }
    double dDydy() const { return dy1; }
    double dDzdz() const { return dz1; }
};

#endif                          // __t3_Test_Prob_hh__

//---------------------------------------------------------------------------//
//                              end of t3/Test_Prob.hh
//---------------------------------------------------------------------------//
