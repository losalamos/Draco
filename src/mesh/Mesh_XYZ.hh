//----------------------------------*-C++-*----------------------------------//
// Mesh_XYZ.hh
// RHS Linux User
// Thu Nov 20 22:08:34 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __mesh_Mesh_XYZ_hh__
#define __mesh_Mesh_XYZ_hh__

#include "Mesh_DB.hh"
#include "ds++/Mat.hh"
#include "c4/NodeInfo.hh"

//===========================================================================//
// class Mesh_XYZ - 

// 
//===========================================================================//

class Mesh_XYZ : private Mesh_DB, private C4::NodeInfo
{
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

  public:
    Mesh_XYZ( const Mesh_DB& mdb );
//     Mesh_XYZ( const Mesh_XYZ& );
//     ~Mesh_XYZ();
//     Mesh_XYZ& operator=( const Mesh_XYZ& );
};

#endif                          // __mesh_Mesh_XYZ_hh__

//---------------------------------------------------------------------------//
//                              end of mesh/Mesh_XYZ.hh
//---------------------------------------------------------------------------//
